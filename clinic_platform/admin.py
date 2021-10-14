from django.contrib import admin

# Register your models here.


from .models import PanelGeneList, CurrentSettings, GeneInfo, Experiment, TranscriptInfo, ExonInfo, CDSInfo, GeneExpression, ExpressionSet, SampleInfo


admin.site.register(PanelGeneList)
admin.site.register(CurrentSettings)
admin.site.register(Experiment)
admin.site.register(GeneInfo)
admin.site.register(TranscriptInfo)
admin.site.register(ExonInfo)
admin.site.register(CDSInfo)
admin.site.register(ExpressionSet)
admin.site.register(GeneExpression)
admin.site.register(SampleInfo)

