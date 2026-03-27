from django import forms
from sympy import sympify, symbols, SympifyError
from .models import SurfaceRecord


class SurfaceRecordForm(forms.ModelForm):
    class Meta:
        model = SurfaceRecord
        fields = ['name', 'X', 'Y', 'Z', 'parameter_names',
                  'u_min', 'u_max', 'v_min', 'v_max', 'u_identify', 'v_identify']

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
        param_str = cleaned.get('parameter_names', '')
        try:
            params = symbols(param_str)
            param_set = set(params) if isinstance(params, tuple) else {params}
        except Exception:
            return cleaned  # parameter_names error already reported

        for field in ('X', 'Y', 'Z'):
            val = cleaned.get(field)
            if not val:
                continue
            try:
                free = sympify(val).free_symbols - param_set
                if free:
                    self.add_error(field, f"Unknown symbol(s): {', '.join(sorted(str(s) for s in free))}")
            except SympifyError:
                pass  # already caught in clean_X/Y/Z

        return cleaned
