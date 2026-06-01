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
                  'domain_type', 'coord_type', 'r_min', 'r_max',
                  'boundary_identify', 'output_type']
        # X/Y/Z are TextField on the model → render as expandable monospace
        # textareas. Defined here (not hardcoded in the template) so the form
        # owns presentation and there is no length cap to keep in sync.
        widgets = {
            f: forms.Textarea(attrs={
                'class': 'form-control surface-expr',
                'rows': 1,
                'style': 'font-family:monospace; resize:vertical; '
                         'overflow:hidden; white-space:pre-wrap; '
                         'word-break:break-all;',
            })
            for f in ('X', 'Y', 'Z')
        }

    def _float_to_display(self, v):
        """Format float for display: use comma decimal, strip trailing zeros."""
        if v is None:
            return ''
        s = f'{v:g}'.replace('.', ',')
        return s

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Rendered as a checkbox (absent = "no"); not required so an unchecked
        # box (which submits nothing) is fine — clean() normalizes it.
        self.fields['boundary_identify'].required = False
        # Pre-fill bound fields with comma-decimal display values from instance
        if self.instance and self.instance.pk:
            for field in ('u_min', 'u_max', 'v_min', 'v_max', 'r_min', 'r_max'):
                val = getattr(self.instance, field, None)
                if val is not None and field not in self.data:
                    self.initial[field] = self._float_to_display(val)
        # The X/Y/Z textarea widgets carry a fixed class from Meta, so add
        # Bootstrap's `is-invalid` on bound errors to surface .invalid-feedback.
        if self.is_bound:
            for name in ('X', 'Y', 'Z'):
                if name in self.errors:
                    w = self.fields[name].widget
                    w.attrs['class'] = (w.attrs.get('class', '') + ' is-invalid').strip()

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
        # Checkbox semantics: value 'an' present = antipodal, absent = none.
        cleaned['boundary_identify'] = (
            'an' if self.data.get('boundary_identify') == 'an' else 'no')
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

        if domain_type == 'rect':
            # Rect domain bounds are (u_min, u_max, v_min, v_max); a zero- or
            # negative-width interval is degenerate and crashes triangulation
            # downstream (a 500). Reject at the form. (Disk ignores u/v bounds —
            # it uses r_min/r_max + 0..2π — so this check is rect-only.)
            u_min = cleaned.get('u_min'); u_max = cleaned.get('u_max')
            v_min = cleaned.get('v_min'); v_max = cleaned.get('v_max')
            if u_min is not None and u_max is not None and u_min >= u_max:
                self.add_error('u_max', "Upper bound must be greater than lower bound (need u_min < u_max).")
            if v_min is not None and v_max is not None and v_min >= v_max:
                self.add_error('v_max', "Upper bound must be greater than lower bound (need v_min < v_max).")

        boundary_identify = cleaned.get('boundary_identify', 'no')

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
            # Antipodal glues the OUTER boundary (u,v)~(-u,-v) → ℝP² (Boy
            # surface). It works at any outer radius R (the involution is
            # σ_R(z) = -R²/z̄), so no r_max constraint; the inner radius is free.
        elif boundary_identify == 'an':
            self.add_error('boundary_identify',
                           "Antipodal boundary identification is only valid for a "
                           "disk domain.")

        return cleaned
