from django.urls import path
from .views import SnakemakeView

urlpatterns = [
    path("run-main-workflow/", SnakemakeView.as_view(), name="run-main-workflow"),
]
