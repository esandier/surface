# On repart sur une base ultra-légère (ubi9-minimal)
FROM image-registry.openshift-image-registry.svc:5000/openshift/python:3.12-ubi9-minimal

USER root
# On installe juste le strict nécessaire
RUN microdnf install -y gcc python3-devel openldap-devel shadow-utils && microdnf clean all

WORKDIR /opt/app-root/src
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .
USER 1001
CMD ["python", "manage.py", "runserver", "0.0.0.0:8080"]