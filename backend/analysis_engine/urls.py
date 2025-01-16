from django.urls import path
from .views import GetUserAnalysisView, CreateUserAnalysisView, ToolUpdateView, ToolVersionCreateView, GetAnalysisDetailView

urlpatterns = [
    path("run-main-workflow/", CreateUserAnalysisView.as_view(), name="run-main-workflow"),
    path("update-tool-version/", ToolVersionCreateView.as_view(), name="update-tool-version"),
    path("run-tool-update-workflow/", ToolUpdateView.as_view(), name="run-tool-update-workflow"),
    path("get-submitted-workflow/", GetUserAnalysisView.as_view(), name="get-submitted-workflow"),
    path("get-specific-workflow/", GetAnalysisDetailView.as_view(), name="get-specific-workflow"),
]
