from django.db import models
from django.urls import reverse

# Create your models here.

class SurfaceRecord(models.Model):
    class Meta:
        ordering = ['name']

    class Identify(models.TextChoices):
        MOEBIUS = 'mo'
        CYLINDRE = 'cy'
        NONE = 'no'
    name = models.CharField(max_length = 100)
    X = models.CharField(max_length=400, default = 'u') 
    Y = models.CharField(max_length=400, default = 'v')
    Z = models.CharField(max_length=400, default = '0')
    parameter_names = models.CharField(max_length=30, default = 'u v')
    u_min = models.FloatField(default = 0)
    u_max = models.FloatField(default = 1)
    v_min = models.FloatField(default = 0)
    v_max = models.FloatField(default = 1)
    u_identify = models.CharField(
        max_length=2,
        choices=Identify.choices,
        default=Identify.NONE,
    )
    v_identify = models.CharField(
        max_length=2,
        choices=Identify.choices,
        default=Identify.NONE,
    )    

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        """Returns the url to access a particular instance of the model."""
        return reverse('surface-play', kwargs={'pk': self.pk})

