from django.contrib import admin
from .models import Transformation, Condition, Reaction

# Register your models here.
admin.site.register(Transformation)
admin.site.register(Condition)
admin.site.register(Reaction)
