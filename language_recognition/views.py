
from django.shortcuts import render
from django.views.generic.base import View
from .forms import DetectLanguageForm
from .processing import vepQuery, wikiPathwayQuery, buildNetwork

class TranslatePhraseView(View):
	template_name = "index.html"

	def get(self, request):
		form = DetectLanguageForm()
		return render(request, self.template_name, {'form': form})

	def post(self, request):
		form = DetectLanguageForm(request.POST)
		countVEP = 0
		if form.is_valid():
			resultVEP, listGene, SNP_GENE = vepQuery(form.data['inputData'])
			if isinstance(SNP_GENE, int):
				return render(request, self.template_name, {'form': form, 
														'resultVEP':resultVEP, 
														'resultWikiPathway':'<h3> </h3>',
														'networkMap': '<h3> </h3>'})
			else:
				resultVEP = '<h3> Queried Genes from VEP</h3>' + resultVEP
				resultWikiPathway, GENE_PATHWAY =  wikiPathwayQuery(listGene)
				resultWikiPathway = '<h3> Queried Pathways from WikiPathway</h3>' + resultWikiPathway

				networkMap = buildNetwork(SNP_GENE, GENE_PATHWAY)
				networkMap = '<h3> SNP-GENE-PATHWAY Network Map</h3>' + networkMap


				return render(request, self.template_name, {'form': form, 
														'resultVEP':resultVEP, 
														'resultWikiPathway':resultWikiPathway,
														'networkMap': networkMap})