#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

namespace {

struct DensityPoint {
  double x;
  double log_density;
  double derivative;
};

struct SamplerCounters {
  double proposals = 0;
  double rejections = 0;
};

double log_phi_plus_exp(double eta, double log_phi) {
  return std::max(eta, log_phi) +
    std::log1p(std::exp(-std::abs(eta - log_phi)));
}

double log_uniform() {
  return std::log(R::runif(0, 1));
}

// Integrate a tangent exp(intercept + slope*x) on its support interval.
double tangent_log_mass(double intercept, double slope, double left, double right) {
  if (!(left < right)) return R_NegInf;
  if (slope == 0) {
    if (!std::isfinite(left) || !std::isfinite(right))
      Rcpp::stop("ARS encountered a non-integrable flat tail.");
    return intercept + std::log(right - left);
  }
  if (slope > 0) {
    if (!std::isfinite(right)) Rcpp::stop("ARS upper tail is not integrable.");
    const double log_difference = std::isfinite(left) ?
      std::log(-std::expm1(slope * (left - right))) : 0;
    return intercept + slope * right + log_difference - std::log(slope);
  }
  if (!std::isfinite(left)) Rcpp::stop("ARS lower tail is not integrable.");
  const double log_difference = std::isfinite(right) ?
    std::log(-std::expm1(slope * (right - left))) : 0;
  return intercept + slope * left + log_difference - std::log(-slope);
}

double sample_tangent(double slope, double left, double right) {
  const double u = R::runif(0, 1);
  if (slope == 0) return left + (right - left) * u;
  if (slope > 0) {
    const double difference = std::isfinite(left) ?
      std::expm1(slope * (left - right)) : -1;
    return right + std::log1p((1 - u) * difference) / slope;
  }
  const double difference = std::isfinite(right) ?
    std::expm1(slope * (right - left)) : -1;
  return left + std::log1p((1 - u) * difference) / slope;
}

// Adaptive rejection sampling uses an upper envelope of supporting tangents.
// A one-sided derivative at beta=0 is a valid supergradient of the concave AL
// log density. The cusp therefore needs no smoothing or Gaussian replacement.
double sample_log_concave(
  const std::function<DensityPoint(double)>& density,
  double initial, double scale, double lower, SamplerCounters& counters
) {
  std::vector<DensityPoint> points;
  double left = std::isfinite(lower) ? lower : std::min(initial, 0.0) - scale;
  double right = std::max(initial, std::isfinite(lower) ? lower : 0.0) + scale;
  DensityPoint left_point = density(left);
  DensityPoint right_point = density(right);
  double expansion = scale;
  for (int attempt = 0; !std::isfinite(lower) && left_point.derivative <= 0; ++attempt) {
    if (attempt >= 100) Rcpp::stop("ARS could not bracket its lower tail.");
    expansion *= 2;
    left -= expansion;
    left_point = density(left);
  }
  expansion = scale;
  for (int attempt = 0; right_point.derivative >= 0; ++attempt) {
    if (attempt >= 100) Rcpp::stop("ARS could not bracket its upper tail.");
    expansion *= 2;
    right += expansion;
    right_point = density(right);
  }
  points.push_back(left_point);
  const double middle = std::max(left, std::min(initial, right));
  if (middle > left && middle < right) points.push_back(density(middle));
  points.push_back(right_point);

  for (int attempt = 0; attempt < 1000; ++attempt) {
    std::vector<DensityPoint> hull;
    std::vector<double> starts;
    for (const DensityPoint& point : points) {
      if (!std::isfinite(point.x) || !std::isfinite(point.log_density) ||
          !std::isfinite(point.derivative)) Rcpp::stop("Non-finite ARS support point.");
      bool redundant = false;
      double start = R_NegInf;
      while (!hull.empty()) {
        const DensityPoint& previous = hull.back();
        const double slope_difference = previous.derivative - point.derivative;
        const double previous_intercept = previous.log_density - previous.derivative * previous.x;
        const double intercept = point.log_density - point.derivative * point.x;
        if (slope_difference < -1e-10 * (1 + std::abs(previous.derivative)))
          Rcpp::stop("ARS conditional log density is not concave.");
        if (slope_difference <= 0) {
          if (intercept >= previous_intercept) {
            redundant = true;
            break;
          }
          hull.pop_back();
          starts.pop_back();
          continue;
        }
        start = (intercept - previous_intercept) / slope_difference;
        if (start > starts.back()) break;
        hull.pop_back();
        starts.pop_back();
      }
      if (!redundant) {
        if (hull.empty()) start = R_NegInf;
        hull.push_back(point);
        starts.push_back(start);
      }
    }
    std::vector<double> masses(hull.size());
    double maximum = R_NegInf;
    for (std::size_t j = 0; j < hull.size(); ++j) {
      const double interval_left = std::max(lower, starts[j]);
      const double interval_right = j + 1 < hull.size() ? starts[j + 1] : R_PosInf;
      const DensityPoint& point = hull[j];
      masses[j] = tangent_log_mass(point.log_density - point.derivative * point.x,
        point.derivative, interval_left, interval_right);
      maximum = std::max(maximum, masses[j]);
    }
    if (!std::isfinite(maximum)) Rcpp::stop("ARS envelope has no finite mass.");
    double total = 0;
    for (double& mass : masses) {
      mass = std::exp(mass - maximum);
      total += mass;
    }
    double threshold = R::runif(0, total);
    std::size_t selected = 0;
    while (selected + 1 < masses.size() && threshold > masses[selected]) {
      threshold -= masses[selected++];
    }
    const double interval_left = std::max(lower, starts[selected]);
    const double interval_right = selected + 1 < hull.size() ? starts[selected + 1] : R_PosInf;
    const DensityPoint& tangent = hull[selected];
    const double value = sample_tangent(tangent.derivative, interval_left, interval_right);
    if (!std::isfinite(value) || value <= lower) continue;
    const DensityPoint candidate = density(value);
    const double upper = tangent.log_density + tangent.derivative * (value - tangent.x);
    const double difference = candidate.log_density - upper;
    if (!std::isfinite(candidate.log_density) || !std::isfinite(candidate.derivative))
      Rcpp::stop("Non-finite ARS proposal density.");
    const double rounding_tolerance = 64 * std::numeric_limits<double>::epsilon() *
      (1 + std::abs(candidate.log_density) + std::abs(upper));
    if (difference > rounding_tolerance)
      Rcpp::stop("ARS supporting envelope fell below the conditional density.");
    counters.proposals += 1;
    if (log_uniform() <= std::min(0.0, difference)) return value;
    counters.rejections += 1;
    const auto position = std::lower_bound(points.begin(), points.end(), value,
      [](const DensityPoint& point, double x) { return point.x < x; });
    if (position == points.end() || position->x != value) points.insert(position, candidate);
  }
  Rcpp::stop("ARS exceeded its rejection limit.");
  return NA_REAL;
}

double positive_normal(double mean, double sd) {
  // Anchor means are nonnegative. Sampling the upper-tail probability avoids
  // rounding an inverse-CDF argument to one for a distant positive draw.
  for (int attempt = 0; attempt < 100; ++attempt) {
    const double tail = R::runif(0, 1) * R::pnorm(mean / sd, 0, 1, true, false);
    const double value = mean + sd * R::qnorm(tail, 0, 1, false, false);
    if (std::isfinite(value) && value > 0) return value;
  }
  Rcpp::stop("Could not draw a finite positive slab auxiliary variable.");
  return NA_REAL;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List tf_activity_gibbs_ars_cpp(
  Rcpp::List common, Rcpp::NumericMatrix eta0, Rcpp::NumericVector phi,
  Rcpp::List settings
) {
  const int cells = Rcpp::as<int>(common["N"]);
  const int genes = Rcpp::as<int>(common["G"]);
  const bool has_condition = Rcpp::as<bool>(common["has_condition"]);
  const Rcpp::NumericMatrix y = Rcpp::as<Rcpp::NumericMatrix>(common["Y"]);
  const Rcpp::NumericVector anchor = common["activity_anchor"];
  const Rcpp::NumericVector shift = common["condition_shift"];
  const Rcpp::NumericVector probability = common["active_prior_probability"];
  const Rcpp::NumericVector lambda_positive = common["beta_lambda_positive"];
  const Rcpp::NumericVector lambda_negative = common["beta_lambda_negative"];
  const Rcpp::NumericVector interaction_sd = common["target_interaction_sd"];
  Rcpp::NumericVector beta = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(common["beta_mean_init"]));
  Rcpp::NumericVector delta = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(common["beta_delta_init"]));
  const double sigma = Rcpp::as<double>(common["sigma_activity"]);
  const double center = Rcpp::as<double>(common["target_center"]);
  const double scale = Rcpp::as<double>(common["target_scale"]);
  const int warmup = Rcpp::as<int>(settings["iter_warmup"]);
  const int samples = Rcpp::as<int>(settings["iter_sampling"]);
  const double move_sd = Rcpp::as<double>(settings["scale_move_sd"]);
  const int iterations = warmup + samples;
  Rcpp::NumericVector auxiliary(cells), activity(cells), standardized(cells), log_phi(genes);
  Rcpp::LogicalVector active(cells);
  for (int g = 0; g < genes; ++g) {
    log_phi[g] = std::log(phi[g]);
    if (!has_condition) delta[g] = 0;
  }
  for (int i = 0; i < cells; ++i) {
    auxiliary[i] = std::max(anchor[i], sigma);
    active[i] = anchor[i] > 0;
    activity[i] = active[i] ? auxiliary[i] : 0;
    standardized[i] = (activity[i] - center) / scale;
  }
  Rcpp::NumericMatrix activity_samples(cells, samples), beta_samples(samples, genes),
    delta_samples(samples, genes);
  Rcpp::NumericVector log_joint_trace(iterations), activity_mean_trace(samples),
    active_count_trace(samples);
  SamplerCounters counters;
  int warmup_moves = 0, retained_moves = 0;

  auto cell_likelihood = [&](int i, double value,
      const Rcpp::NumericVector& beta_value, const Rcpp::NumericVector& delta_value) {
    double result = 0;
    const double predictor = (value - center) / scale;
    for (int g = 0; g < genes; ++g) {
      const double eta = eta0(i, g) + predictor * (beta_value[g] + shift[i] * delta_value[g]);
      result += y(i, g) * eta - (y(i, g) + phi[g]) * log_phi_plus_exp(eta, log_phi[g]);
    }
    return result;
  };
  auto joint_density = [&](const Rcpp::NumericVector& beta_value,
      const Rcpp::NumericVector& delta_value, const Rcpp::NumericVector& activity_value) {
    double result = 0;
    for (int g = 0; g < genes; ++g) {
      result += beta_value[g] >= 0 ? -lambda_positive[g] * beta_value[g] :
        lambda_negative[g] * beta_value[g];
      if (has_condition) result -= 0.5 * std::pow(delta_value[g] / interaction_sd[g], 2);
    }
    for (int i = 0; i < cells; ++i) {
      result += cell_likelihood(i, activity_value[i], beta_value, delta_value);
      // The inactive auxiliary slab value is part of the augmented target.
      // Normalizing constants are fixed within this nuisance component.
      const double slab_value = active[i] ? activity_value[i] : auxiliary[i];
      result += active[i] ? std::log(probability[i]) : std::log1p(-probability[i]);
      result -= 0.5 * std::pow((slab_value - anchor[i]) / sigma, 2);
    }
    return result;
  };

  for (int iteration = 0; iteration < iterations; ++iteration) {
    Rcpp::checkUserInterrupt();
    for (int g = 0; g < genes; ++g) {
      auto coefficient_density = [&](double value, bool interaction) {
        double log_density = interaction ? -0.5 * std::pow(value / interaction_sd[g], 2) :
          (value >= 0 ? -lambda_positive[g] * value : lambda_negative[g] * value);
        double derivative = interaction ? -value / std::pow(interaction_sd[g], 2) :
          (value >= 0 ? -lambda_positive[g] : lambda_negative[g]);
        for (int i = 0; i < cells; ++i) {
          const double beta_value = interaction ? beta[g] : value;
          const double delta_value = interaction ? value : delta[g];
          const double eta = eta0(i, g) + standardized[i] * (beta_value + shift[i] * delta_value);
          const double predictor = standardized[i] * (interaction ? shift[i] : 1);
          log_density += y(i, g) * eta - (y(i, g) + phi[g]) * log_phi_plus_exp(eta, log_phi[g]);
          derivative += predictor * (y(i, g) - (y(i, g) + phi[g]) *
            R::plogis(eta - log_phi[g], 0, 1, true, false));
        }
        return DensityPoint{value, log_density, derivative};
      };
      const double beta_scale = std::hypot(1 / lambda_positive[g], 1 / lambda_negative[g]);
      beta[g] = sample_log_concave([&](double value) {
        return coefficient_density(value, false);
      }, beta[g], beta_scale, R_NegInf, counters);
      if (has_condition) {
        delta[g] = sample_log_concave([&](double value) {
          return coefficient_density(value, true);
        }, delta[g], interaction_sd[g], R_NegInf, counters);
      }
    }
    for (int i = 0; i < cells; ++i) {
      if (active[i]) {
        auto activity_density = [&](double value) {
          double log_density = -0.5 * std::pow((value - anchor[i]) / sigma, 2);
          double derivative = -(value - anchor[i]) / (sigma * sigma);
          for (int g = 0; g < genes; ++g) {
            const double coefficient = (beta[g] + shift[i] * delta[g]) / scale;
            const double eta = eta0(i, g) + coefficient * (value - center);
            log_density += y(i, g) * eta - (y(i, g) + phi[g]) * log_phi_plus_exp(eta, log_phi[g]);
            derivative += coefficient * (y(i, g) - (y(i, g) + phi[g]) *
              R::plogis(eta - log_phi[g], 0, 1, true, false));
          }
          return DensityPoint{value, log_density, derivative};
        };
        auxiliary[i] = sample_log_concave(activity_density, auxiliary[i], sigma, 0, counters);
      } else {
        auxiliary[i] = positive_normal(anchor[i], sigma);
      }
      const double log_odds = std::log(probability[i]) - std::log1p(-probability[i]) +
        cell_likelihood(i, auxiliary[i], beta, delta) - cell_likelihood(i, 0, beta, delta);
      active[i] = R::runif(0, 1) < R::plogis(log_odds, 0, 1, true, false);
      activity[i] = active[i] ? auxiliary[i] : 0;
      standardized[i] = (activity[i] - center) / scale;
    }

    // Symmetric log-scale proposal on the augmented state. The affine activity
    // transformation preserves the centered beta*A predictor for active cells.
    // Inactive cells stay at A=0, and their changed likelihood is evaluated too.
    const double log_multiplier = R::rnorm(0, move_sd);
    const double multiplier = std::exp(log_multiplier);
    Rcpp::NumericVector proposed_beta = beta * multiplier;
    Rcpp::NumericVector proposed_delta = delta * multiplier;
    Rcpp::NumericVector proposed_activity = Rcpp::clone(activity);
    int active_count = 0;
    bool valid_move = std::isfinite(multiplier) && multiplier > 0;
    for (int i = 0; i < cells; ++i) {
      if (!active[i]) continue;
      ++active_count;
      proposed_activity[i] = center + (activity[i] - center) / multiplier;
      valid_move = valid_move && std::isfinite(proposed_activity[i]) && proposed_activity[i] > 0;
    }
    const double current_joint = joint_density(beta, delta, activity);
    double final_joint = current_joint;
    if (valid_move) {
      const double proposed_joint = joint_density(proposed_beta, proposed_delta, proposed_activity);
      const int coefficient_count = genes * (has_condition ? 2 : 1);
      const double log_acceptance = proposed_joint - current_joint +
        (coefficient_count - active_count) * log_multiplier;
      if (std::isfinite(proposed_joint) && log_uniform() < std::min(0.0, log_acceptance)) {
        beta = proposed_beta;
        delta = proposed_delta;
        activity = proposed_activity;
        final_joint = proposed_joint;
        for (int i = 0; i < cells; ++i) {
          if (active[i]) auxiliary[i] = activity[i];
          standardized[i] = (activity[i] - center) / scale;
        }
        if (iteration < warmup) ++warmup_moves; else ++retained_moves;
      }
    }
    log_joint_trace[iteration] = final_joint;
    if (iteration >= warmup) {
      const int draw = iteration - warmup;
      double activity_total = 0;
      for (int i = 0; i < cells; ++i) {
        activity_samples(i, draw) = activity[i];
        activity_total += activity[i];
      }
      for (int g = 0; g < genes; ++g) {
        beta_samples(draw, g) = beta[g];
        delta_samples(draw, g) = delta[g];
      }
      activity_mean_trace[draw] = activity_total / cells;
      active_count_trace[draw] = active_count;
    }
  }
  return Rcpp::List::create(
    Rcpp::_["activity_samples"] = activity_samples,
    Rcpp::_["beta_samples"] = beta_samples,
    Rcpp::_["delta_samples"] = delta_samples,
    Rcpp::_["log_joint_trace"] = log_joint_trace,
    Rcpp::_["activity_mean_trace"] = activity_mean_trace,
    Rcpp::_["active_count_trace"] = active_count_trace,
    Rcpp::_["ars_proposals"] = counters.proposals,
    Rcpp::_["ars_rejections"] = counters.rejections,
    Rcpp::_["scale_accepted_warmup"] = warmup_moves,
    Rcpp::_["scale_accepted_sampling"] = retained_moves
  );
}
