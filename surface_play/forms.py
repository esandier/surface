from django import forms
from sympy import sympify, symbols, SympifyError, pi, E, N as sympyN
from .models import SurfaceRecord


class SympyFloatField(forms.CharField):
    """Text field that accepts SymPy expressions (e.g. '2*pi') and returns a float.
    Accepts comma as decimal separator (French convention)."""

    def to_python(self, value):
        value = super().to_python(value)
        if not value:
            return None
        # Accept comma as decimal separator
        value = value.replace(',', '.')
        try:
            result = sympify(value, locals={'pi': pi, 'E': E})
            f = float(sympyN(result))
        except Exception:
            raise forms.ValidationError("Enter a number or a mathematical expression (e.g. 2*pi).")
        if not (f == f):  # NaN check
            raise forms.ValidationError("Expression must evaluate to a finite number.")
        return f


class SurfaceRecordForm(forms.ModelForm):

    u_min = SympyFloatField(initial='0')
    u_max = SympyFloatField(initial='1')
    v_min = SympyFloatField(initial='0')
    v_max = SympyFloatField(initial='1')
    r_min = SympyFloatField(initial='0')
    r_max = SympyFloatField(initial='1')

    class Meta:
        model = SurfaceRecord
        fields = ['name', 'X', 'Y', 'Z', 'parameter_names',
                  'u_min', 'u_max', 'v_min', 'v_max', 'u_identify', 'v_identify',
                  'domain_type', 'coord_type', 'r_min', 'r_max', 'output_type']

    def _float_to_display(self, v):
        """Format float for display: use comma decimal, strip trailing zeros."""
        if v is None:
            return ''
        s = f'{v:g}'.replace('.', ',')
        return s

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Pre-fill bound fields with comma-decimal display values from instance
        if self.instance and self.instance.pk:
            for field in ('u_min', 'u_max', 'v_min', 'v_max', 'r_min', 'r_max'):
                val = getattr(self.instance, field, None)
                if val is not None and field not in self.data:
                    self.initial[field] = self._float_to_display(val)

    def _check_expr(self, field_name):
        val = self.cleaned_data.get(field_name, '')
        try:
            sympify(val)
        except SympifyError as e:
            raise forms.ValidationError(f"Invalid expression: {e}")
        return val

    def clean_X(self): return self._check_expr('X')
    def clean_Y(self): return self._check_expr('Y')
    def clean_Z(self): return self._check_expr('Z')

    def clean(self):
        cleaned = super().clean()
        domain_type = cleaned.get('domain_type', 'rect')
        coord_type  = cleaned.get('coord_type', 'ca')
        param_str   = cleaned.get('parameter_names', '')

        try:
            params = symbols(param_str)
            param_list = list(params) if isinstance(params, tuple) else [params]
        except Exception:
            return cleaned

        if len(param_list) != 2:
            self.add_error('parameter_names',
                           f"Expected 2 variable names, got {len(param_list)}.")
            return cleaned

        param_set = set(param_list)

        if domain_type == 'disk' and coord_type == 'po':
            return cleaned

        for field in ('X', 'Y', 'Z'):
            val = cleaned.get(field)
            if not val:
                continue
            try:
                free = sympify(val).free_symbols - param_set
                if free:
                    self.add_error(field,
                                   f"Unknown symbol(s): {', '.join(sorted(str(s) for s in free))}")
            except SympifyError:
                pass

        if domain_type == 'disk':
            r_min = cleaned.get('r_min', 0.0)
            r_max = cleaned.get('r_max', 1.0)
            if r_min is not None and r_max is not None:
                if r_min < 0:
                    self.add_error('r_min', "Inner radius must be ≥ 0.")
                if r_max <= 0:
                    self.add_error('r_max', "Outer radius must be > 0.")
                if r_min >= r_max:
                    self.add_error('r_max', "Outer radius must be greater than inner radius.")

        return cleaned
