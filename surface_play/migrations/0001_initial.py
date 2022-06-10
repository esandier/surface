# Generated by Django 4.0.4 on 2022-05-26 21:01

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='SurfaceRecord',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=100)),
                ('X', models.CharField(max_length=400)),
                ('Y', models.CharField(max_length=400)),
                ('Z', models.CharField(max_length=400)),
                ('paramNames', models.CharField(max_length=30)),
                ('u_identify', models.CharField(choices=[('mo', 'Moebius'), ('cy', 'Cylindre'), ('no', 'None')], default='no', max_length=2)),
                ('v_identify', models.CharField(choices=[('mo', 'Moebius'), ('cy', 'Cylindre'), ('no', 'None')], default='no', max_length=2)),
            ],
        ),
    ]