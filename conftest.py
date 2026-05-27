"""Repo-level pytest config: bootstrap Django so test modules can import
ORM models / use the test client without each having to call
`django.setup()` themselves.

Only takes effect when running under pytest. Leaves non-Django tests
(Layer C/O units, pipeline) unaffected — they just never touch the ORM.
"""

import os

import django


os.environ.setdefault("DJANGO_SETTINGS_MODULE", "surfaces_project.settings")
django.setup()
