# -*- coding: utf-8 -*-

import sys


def runme(count_matrix_path, design_matrix_path, gene_column):
    import pandas as pd
    import rpy2
    from rpy2.robjects import pandas2ri, Formula, r

    assert rpy2.__version__, '2.9.1' ' Please install rpy2 2.9.1 to run this script'
    assert pd.__version__, '0.19' ' Please install pandas 0.19 to run this script'

    pandas2ri.activate()
    from rpy2.robjects.packages import importr

    try:
        deseq = importr('DESeq2')
    except:
        EnvironmentError('Please install DESeq2 in your R environment')

    # Necessary to translate R dataframe back to Pandas
    to_dataframe = r('function(x) data.frame(x)')

    print('Loading data with pandas')

    count_matrix_df = pd.read_csv(
        count_matrix_path,
        sep=',',
    )

    design_matrix_df = pd.read_csv(
        design_matrix_path,
        sep=',',
        index_col=0
    )

    class py_DESeq2:

        def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
            try:
                assert gene_column in count_matrix.columns, 'Wrong gene id column name'
                gene_id = count_matrix[gene_column]
            except AttributeError:
                sys.exit('Wrong Pandas dataframe?')

            self.dds = None
            self.deseq_result = None
            self.resLFC = None
            self.comparison = None
            self.normalized_count_matrix = None
            self.gene_column = gene_column
            self.gene_id = count_matrix[self.gene_column]

            count_matrix = count_matrix.drop(gene_column, axis=1)

            print(
                f'Number of columns in counts data {count_matrix.shape[1]} | '
                f'Number of rows in design matrix {design_matrix.shape[0]}'
            )

            # Load dataframe into R environment
            # Important: Change to r.data() if you use numpys and rpy2 latests versions
            count_matrix = pandas2ri.py2ri(count_matrix)

            # Assign columns to NULL
            count_matrix.names = rpy2.rinterface.NULL

            self.count_matrix = count_matrix

            self.design_matrix = pandas2ri.py2ri(design_matrix)

            self.design_formula = Formula(design_formula)

        def run_deseq(self, **kwargs):
            self.dds = deseq.DESeqDataSetFromMatrix(
                countData=self.count_matrix,
                colData=self.design_matrix,
                design=self.design_formula
            )
            self.dds = deseq.DESeq(self.dds, **kwargs)
            # Previous script had "deseq.counts" instead
            self.normalized_count_matrix = deseq.counts_DESeqDataSet(self.dds, normalized=True)

        def get_deseq_result(self, **kwargs):

            self.comparison = deseq.resultsNames(self.dds)

            self.deseq_result = deseq.results(self.dds, **kwargs)
            self.deseq_result = to_dataframe(self.deseq_result)
            self.deseq_result = pandas2ri.ri2py(self.deseq_result)  ## back to pandas dataframe
            self.deseq_result[self.gene_column] = self.gene_id.values
            return self.deseq_result

    print('Creating R objects')
    deseq2_exp = py_DESeq2(
        count_matrix=count_matrix_df,
        design_matrix=design_matrix_df,
        design_formula='~ class_label',
        gene_column=gene_column
    )

    print('Running DESeq2 scripts...please be patient')

    deseq2_exp.run_deseq()

    print('Almost done...getting the results ready')

    results = deseq2_exp.get_deseq_result()

    results.to_csv('results.csv')

    print('Done!')


if __name__ == "__main__":
    # count_matrix. "lihc_read_count_duplicates_removed.txt"
    # design_matrix "LIHC_transposed_labels.txt"
    runme(
        '/home/ddomingofernandez/Projects/pathwayforte-resources/data/rnaseq/lihc/lihc_read_count_duplicates_removed.txt',
        '/home/ddomingofernandez/Projects/pathwayforte-resources/data/rnaseq/lihc/LIHC_transposed_labels.txt',
        'gene_symbol'
    )
