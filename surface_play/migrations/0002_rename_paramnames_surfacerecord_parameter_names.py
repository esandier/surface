# Generated by Django 4.0.4 on 2022-05-26 21:55

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('surface_play', '0001_initial'),
    ]

    operations = [
        migrations.RenameField(
            model_name='surfacerecord',
            old_name='paramNames',
            new_name='parameter_names',
        ),
    ]
