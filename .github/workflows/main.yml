name: Project 2 Python Tests

on:
  push:
  pull_request:
  workflow_dispatch:

jobs:
  test:
    runs-on: ubuntu-latest

    defaults:
      run:
        working-directory: ammonia_absorption_chiller

    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.11"

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          python -m pip install -r requirements.txt
          python -m pip install pytest
      - name: Verify NumPy
        run: python -c "import numpy; print(numpy.__version__)"

      - name: Run tests
        run: python -m pytest -v
