"""Optional, evidence-grounded language-model explanations."""

from .llm_explainer import explain_predictions, llm_available

__all__ = ["explain_predictions", "llm_available"]
