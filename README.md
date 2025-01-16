# XPLORECoV

XPLORECoV is your gateway to advanced SARS-CoV-2 genomic analysis. Combining cutting-edge bioinformatics workflows with AI-driven intelligence, this free, no-login-required platform transforms raw genomic data into actionable epidemiological insights. Say goodbye to outdated filters! XPLORECoV’s powerful AI capabilities handle complex, multifactorial questions with ease, delivering intuitive results on mutation patterns, lineage trends, and spatiotemporal epidemiology. Whether you're an academic, a researcher, or a public health expert, XPLORECoV empowers you to dive deeper with high-throughput pipelines, global datasets, and detailed visual reports

# XploreCoV Setup Guide

This guide provides step-by-step instructions to set up both the backend and frontend of the XploreCoV application.

## Backend Setup

Follow these steps to set up the backend of the project:

### Prerequisites

-   **Python 3.11** must be installed on your system.
-   Install **Poetry** for dependency management.
-   Install **Micromamba** for workflow tools.

### Steps

1. **Install Dependencies**:

    - Navigate to the backend directory:
        ```bash
        cd backend
        ```
    - Install dependencies using Poetry:
        ```bash
        poetry install
        ```
    - Activate the Poetry environment:
        ```bash
        poetry shell
        ```

2. **Create a `.env` File**:

    - Inside the `backend` folder, create a file named `.env` with the following keys:
        ```env
        BASE_HOST=<your_base_host>
        BASE_PORT=<your_base_port>
        DEBUG=<true_or_false>
        BASE_URL=<your_base_url>
        BACKEND_WEBSOCKET_UUID=<uuid_placeholder>
        ```

3. **Apply Migrations**:

    - Run the following commands to set up the database:
        ```bash
        python manage.py makemigrations
        python manage.py migrate
        ```

4. **Create a Superuser**:

    - Create a superuser account:
        ```bash
        python manage.py createsuperuser
        ```

5. **Generate Backend UUID**:

    - Generate a UUID for the backend:
        ```bash
        python manage.py create_backend_uuid
        ```
    - Add the generated UUID to the `.env` file under `BACKEND_WEBSOCKET_UUID`.

6. **Install Workflow Tools**:

    - Run the following command to install the workflow tools:
        ```bash
        micromamba install -r .workflow-venv -n xplorecov
        ```

7. **Run Required Services**:

    - Start Redis:
        ```bash
        redis-server
        ```
        Ensure that logs are stored in the `backend/database` folder.
    - Start Celery:
        ```bash
        celery -A xplorecov_backend worker -l info -E
        ```
    - Start Gunicorn:
        ```bash
        gunicorn xplorecov_backend.asgi:application -c gunicorn.config.py
        ```

8. **Update Workflow Tools**:

    - Update tools and fetch resources by calling the API:
        ```
        xplorecov/api/job/run-tool-update-workflow/
        ```

9. **Folder Requirements**:
    - Ensure the following folders exist:
        - `datalake`
        - `database`
    - If these folders do not exist, create them before running any commands.

---

## Frontend Setup

Follow these steps to set up the frontend of the project:

### Steps

1. **Install Dependencies**:

    - Navigate to the frontend directory:
        ```bash
        cd frontend
        ```
    - Install dependencies using Bun:
        ```bash
        bun install
        ```

2. **Create a `.env` File**:
    - Inside the `frontend` folder, create a file named `.env` with the following keys:
        ```env
        API_BASE_URL=<api_base_url>
        ROUTER_BASE=<router_base>
        ```

---

## Summary

-   Set up the backend by installing dependencies, creating a `.env` file, applying migrations, and running services.
-   Generate a backend UUID and include it in the `.env` file.
-   Start Redis, Celery, and Gunicorn services.
-   Set up the frontend by installing dependencies and configuring the `.env` file.

With everything set up, the application should be ready to run. Ensure all services and workflows are properly configured before starting the application.
