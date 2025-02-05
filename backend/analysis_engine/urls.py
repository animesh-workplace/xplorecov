from django.urls import path
from .views import (
    ToolUpdateView,
    AddBulkReportsView,
    GetUserAnalysisView,
    AddChatMessagesView,
    ToolVersionCreateView,
    GetAnalysisDetailView,
    CreateUserAnalysisView,
)

urlpatterns = [
    path(
        "run-main-workflow/", CreateUserAnalysisView.as_view(), name="run-main-workflow"
    ),
    path(
        "update-tool-version/",
        ToolVersionCreateView.as_view(),
        name="update-tool-version",
    ),
    path(
        "run-tool-update-workflow/",
        ToolUpdateView.as_view(),
        name="run-tool-update-workflow",
    ),
    path(
        "get-submitted-workflow/",
        GetUserAnalysisView.as_view(),
        name="get-submitted-workflow",
    ),
    path(
        "get-specific-workflow/",
        GetAnalysisDetailView.as_view(),
        name="get-specific-workflow",
    ),
    path("reports/create/", AddBulkReportsView.as_view(), name="bulk-create-reports"),
    path(
        "messages/create/",
        AddChatMessagesView.as_view(),
        name="bulk-create-messages",
    ),
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
