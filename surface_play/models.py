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

    class BoundaryIdentify(models.TextChoices):
        ANTIPODAL = 'an'
        NONE = 'no'

    class DomainType(models.TextChoices):
        RECTANGLE = 'rect'
        DISK = 'disk'

    class CoordType(models.TextChoices):
        CARTESIAN = 'ca'
        POLAR = 'po'

    class OutputType(models.TextChoices):
        CARTESIAN = 'ca'
        CYLINDRICAL = 'cy'

    name = models.CharField(max_length = 100)
    # TextField (no length cap) — expressions can be long (e.g. Boy surface).
    # The form renders these as expandable textareas (see SurfaceRecordForm).
    X = models.TextField(default='u')
    Y = models.TextField(default='v')
    Z = models.TextField(default='0')
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
    domain_type = models.CharField(
        max_length=4,
        choices=DomainType.choices,
        default=DomainType.RECTANGLE,
    )
    # Disk-only: 'an' glues the boundary circle antipodally ((u,v)~(-u,-v)),
    # making the disk an ℝP² immersion (e.g. Boy surface). Ignored for rect.
    boundary_identify = models.CharField(
        max_length=2,
        choices=BoundaryIdentify.choices,
        default=BoundaryIdentify.NONE,
    )
    coord_type = models.CharField(
        max_length=2,
        choices=CoordType.choices,
        default=CoordType.CARTESIAN,
    )
    r_min = models.FloatField(default=0.0)
    r_max = models.FloatField(default=1.0)
    output_type = models.CharField(
        max_length=2,
        choices=OutputType.choices,
        default=OutputType.CARTESIAN,
    )
    thumbnail = models.TextField(blank=True, default='')
    initial_elev   = models.FloatField(default=35.0)
    initial_azim   = models.FloatField(default=45.0)
    initial_inplane = models.FloatField(default=0.0)
    initial_zoom   = models.FloatField(default=1.0)
    initial_perspective = models.BooleanField(default=False)

    def __str__(self):
        return self.name

    def get_absolute_url(self):
        """Returns the url to access a particular instance of the model."""
        return reverse('surface-play', kwargs={'pk': self.pk})

