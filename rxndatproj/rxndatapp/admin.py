from django.contrib import admin
from .models import Transformation, Chemical, Condition, Substrate, Reaction

# Register your models here.
admin.site.register(Transformation)
admin.site.register(Chemical)
admin.site.register(Condition)
admin.site.register(Substrate)
admin.site.register(Reaction)