FROM image-registry.openshift-image-registry.svc:5000/openshift/python:3.11-ubi8

WORKDIR /opt/app-root/src

# Copy requirements first so the pip layer is cached on code-only changes
COPY requirements.txt .
RUN chown 1001:0 requirements.txt

USER 1001
RUN pip install --no-cache-dir -r requirements.txt

USER root
COPY . .
RUN chown -R 1001:0 /opt/app-root/src

USER 1001
RUN python manage.py collectstatic --noinput

CMD ["gunicorn", "--bind", "0.0.0.0:8080", "--workers", "2", "surfaces_project.wsgi:application"]
