from django.urls import path
from .views import UserAnalysisView, ToolUpdateView, ToolVersionCreateView

urlpatterns = [
    path("run-main-workflow/", UserAnalysisView.as_view(), name="run-main-workflow"),
    path("update-tool-version/", ToolVersionCreateView.as_view(), name="update-tool-version"),
    path("run-tool-update-workflow/", ToolUpdateView.as_view(), name="run-tool-update-workflow"),
]
