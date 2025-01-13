from django.urls import path
from .views import SnakemakeView, UserAnalysisView

urlpatterns = [
    path("run-main-workflow/", UserAnalysisView.as_view(), name="run-main-workflow"),
]
