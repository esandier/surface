from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('surface_play', '0007_surfacerecord_initial_view'),
    ]

    operations = [
        migrations.AddField(
            model_name='surfacerecord',
            name='initial_perspective',
            field=models.BooleanField(default=False),
        ),
    ]
