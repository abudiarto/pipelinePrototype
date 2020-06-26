from django import forms

class DetectLanguageForm(forms.Form):
    inputData = forms.CharField(
                            required=True,
                            widget=forms.Textarea,
                            error_messages={'required': 'Message is required'}
                            )


    inputData.widget.attrs['class'] = 'form-control'