import scConnect as cn
import scanpy as sc

adata = sc.datasets.paul15()
sc.pp.log1p(adata)


def test_import():
    assert 1+1 == 2

def test_gene_call(adata):
    cn.genecall.meanExpression(adata, groupby="paul15_clusters", transformation="log1p")
    assert list(adata.uns["gene_call"].keys()) == ['1Ery', '2Ery', '3Ery', '4Ery', '5Ery', '6Ery', '7MEP', '8Mk', '9GMP', '10GMP', '11DC', '12Baso', '13Baso', '14Mo', '15Mo', '16Neu', '17Neu', '18Eos', '19Lymph']