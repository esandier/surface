#!/bin/bash
# Install dependencies
pip3 install -r requirements.txt

# Prepare static files
python3 manage.py collectstatic --noinput

# Start the server
python3 -m gunicorn surfaces_project.wsgi:application --bind 0.0.0.0:8000