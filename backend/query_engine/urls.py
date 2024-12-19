from django.urls import path
from .views import VirusDataUploadView, AIModelTrainView, AIModelQueryView

urlpatterns = [
    path('query-ai-model/', AIModelQueryView.as_view(), name='query-ai-model'),
    path('train-ai-model/', AIModelTrainView.as_view(), name='train-ai-model'),
    path('upload-virus-data/', VirusDataUploadView.as_view(), name='upload-virus-data'),
]