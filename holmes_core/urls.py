"""holmes_core URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.conf.urls import include, url
from django.conf.urls.static import static
from django.conf import settings

from clinic_platform import views

urlpatterns = [
    # Vanilla
    url(r'^$', views.home, name='home'),

    url(r'^uploads/$', views.uploads, name='uploads'),
    url(r'^viewGene/(?P<gene_id>\w{0,50})/$', views.gene_view, name='geneView'),
    url(r'^searchGenes/$', views.gene_search, name='geneSearch'),
    url(r'^subNetworks/$', views.sub_graphs, name='subNetworks'),

    # Variant views
    url(r'^varOverview/', views.variant_overview, name='varOverview'),
    url(r'^varOverview/(?P<variant>\w{0,50})/$', views.variant_overview, name='varOverview'),

    # Add panel views
    url(r'^addPanelVCF/$', views.add_vcf_to_panel, name='addPanelVCF'),

    url(r'^subNetworks/Subnet_*?', views.sub_graph_detail, name='subNetworkDetail'),
    url(r'^coverageSummaryGene/$', views.coverage_summary_gene, name='coverageSummaryGene'),
    url(r'^coverageSummarySample/$', views.coverage_summary_sample, name='coverageSummarySample'),
    url(r'^viewRegion/$', views.view_region, name='viewRegion'),
    url(r'^getGeneSubset/$', views.create_subset_file, name='getGeneSubset'),
    url(r'^inputGeneList/$', views.gene_list_input, name='geneListInput'),
    url(r'^genelist2db/$', views.gene_list_to_db, name='genelist2db'),
    url('admin/', admin.site.urls),
    url('process_bam/$', views.process_sample_coverage, name='processSampleCoverage'),
    url('custom_lists/$', views.gene_file_to_db, name='geneFileToDB'),
    url('panelDetails/$', views.panel_detail_view, name='panelDetails'),
    url('addExperiment/$', views.add_experiment, name='addExperiment'),
    url('addExpressionData/$', views.add_expression_data, name='addExpressionData'),

    # Help pages
    url('helpGeneInput/$', views.help_gene_input, name='helpGeneInput'),
    url('helpVCFInput/$', views.help_vcf_input, name='helpVCFInput'),
    url('helpPanelInput/$', views.help_panel_input, name='helpPanelInput'),

] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

#url(r'^(varOverview/?P<variant>\w+)/$', views.variant_overview, name='varOverview'),
#url(r'^varOverview/$', views.variant_overview, name='default_varOverview'),