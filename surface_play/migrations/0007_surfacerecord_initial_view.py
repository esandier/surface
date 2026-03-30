from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('surface_play', '0006_surfacerecord_domain'),
    ]

    operations = [
        migrations.AddField(
            model_name='surfacerecord',
            name='output_type',
            field=models.CharField(
                choices=[('ca', 'Cartesian'), ('cy', 'Cylindrical')],
                default='ca',
                max_length=2,
            ),
        ),
        migrations.AddField(
            model_name='surfacerecord',
            name='initial_elev',
            field=models.FloatField(default=35.0),
        ),
        migrations.AddField(
            model_name='surfacerecord',
            name='initial_azim',
            field=models.FloatField(default=45.0),
        ),
        migrations.AddField(
            model_name='surfacerecord',
            name='initial_inplane',
            field=models.FloatField(default=0.0),
        ),
        migrations.AddField(
            model_name='surfacerecord',
            name='initial_zoom',
            field=models.FloatField(default=1.0),
        ),
    ]
