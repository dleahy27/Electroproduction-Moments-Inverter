#include "Detail.h"

#include "Math/Factory.h"

#include <cmath>
#include <stdexcept>

namespace emi::detail {

MinimizerResult Minimize(const std::shared_ptr<EvaluationContext>& context,
                         const std::vector<double>& start,
                         bool runHesse) {
  if (!context) throw std::invalid_argument("Minimize received a null context");
  if (start.size() != context->freeToFull.size()) {
    throw std::invalid_argument("Start point has the wrong dimension");
  }

  auto analyticObjective = context->cfg.useNumericalGradients
                               ? nullptr
                               : MakeAnalyticChi2(context);
  auto numericalObjective = context->cfg.useNumericalGradients
                                ? MakeNumericalChi2(context)
                                : nullptr;
  std::unique_ptr<ROOT::Math::Minimizer> minimizer(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  if (!minimizer) throw std::runtime_error("Could not create Minuit2/Migrad");

  minimizer->SetMaxFunctionCalls(context->cfg.maxCalls);
  minimizer->SetMaxIterations(context->cfg.maxIterations);
  minimizer->SetTolerance(context->cfg.tolerance);
  minimizer->SetStrategy(context->cfg.strategy);
  minimizer->SetPrintLevel(context->cfg.printLevel);
  minimizer->SetErrorDef(1.0);
  if (numericalObjective) {
    minimizer->SetFunction(*numericalObjective);
  } else {
    minimizer->SetFunction(*analyticObjective);
  }

  for (unsigned i = 0; i < start.size(); ++i) {
    const auto& parameter = context->fullPars[context->freeToFull[i]];
    minimizer->SetLimitedVariable(i, parameter.name.c_str(), start[i],
                                  parameter.step, parameter.low, parameter.high);
  }

  MinimizerResult result;
  result.fitOk = minimizer->Minimize();
  result.hesseOk = runHesse && result.fitOk ? minimizer->Hesse() : false;
  result.status = minimizer->Status();
  result.covarianceStatus = minimizer->CovMatrixStatus();
  result.edm = minimizer->Edm();
  result.chi2 = minimizer->MinValue();
  result.valid = std::isfinite(result.chi2) && result.chi2 < 1e299;
  result.freeValues.assign(minimizer->X(), minimizer->X() + start.size());
  result.minimizer = std::move(minimizer);
  return result;
}

} // namespace emi::detail
