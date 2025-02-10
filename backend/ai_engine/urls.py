from django.urls import path
from .views import (
    AddChatMessagesView,
)

urlpatterns = [
    path("ask-ai/", AddChatMessagesView.as_view(), name="create-messages"),
]

# urlpatterns = [
#     # User Analysis APIs
#     path(
#         "user-analysis/create/",
#         CreateUserAnalysisView.as_view(),
#         name="create-user-analysis",
#     ),
#     path(
#         "user-analysis/list/", GetUserAnalysisView.as_view(), name="list-user-analyses"
#     ),
#     path(
#         "user-analysis/detail/",
#         GetAnalysisDetailView.as_view(),
#         name="detail-user-analysis",
#     ),
#     # Tool Version APIs
#     path(
#         "tool-version/update/",
#         ToolVersionCreateView.as_view(),
#         name="update-tool-version",
#     ),
#     path(
#         "tool/update-workflow/", ToolUpdateView.as_view(), name="update-tool-workflow"
#     ),
#     # Report APIs
#     path(
#         "reports/bulk-create/", AddBulkReportsView.as_view(), name="bulk-create-reports"
#     ),
#     # Chat Message APIs
#     path(
#         "messages/bulk-create/",
#         AddChatMessagesView.as_view(),
#         name="bulk-create-messages",
#     ),
# ]
