from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('surface_play', '0005_surfacerecord_thumbnail'),
    ]

    operations = [
        migrations.AddField(
            model_name='surfacerecord',
            name='domain_type',
            field=models.CharField(
                choices=[('rect', 'Rectangle'), ('disk', 'Disk')],
                default='rect',
                max_length=4,
            ),
        ),
        migrations.AddField(
            model_name='surfacerecord',
            name='coord_type',
            field=models.CharField(
                choices=[('ca', 'Cartesian'), ('po', 'Polar')],
                default='ca',
                max_length=2,
            ),
        ),
        migrations.AddField(
            model_name='surfacerecord',
            name='r_min',
            field=models.FloatField(default=0.0),
        ),
        migrations.AddField(
            model_name='surfacerecord',
            name='r_max',
            field=models.FloatField(default=1.0),
        ),
    ]
