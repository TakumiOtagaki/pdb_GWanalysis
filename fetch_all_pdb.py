from Bio.PDB.PDBList import PDBList

def load_pdb_ids():
    # RNAを含むPDB IDのリストを取得（仮の例として、ここでは具体的なクエリ方法を示しますが、
    # 実際にはPDBのAPIや他のデータベースを利用してリストを取得する必要があります）
    # ここでは、フェッチするPDB IDを手動で指定します。
    ret = []
    with open("/work/gs58/s58007/rna_GWclustering/all_pdb_id.csv") as f:
        for line in f:
            ret.append(line.strip())
    return ret

def fetch_all_rna_structures():
    # PDBListオブジェクトを作成
    pdbl = PDBList()
    
    # RNAを含むPDB IDのリストを取得（仮の例として、ここでは具体的なクエリ方法を示しますが、
    # 実際にはPDBのAPIや他のデータベースを利用してリストを取得する必要があります）
    # ここでは、フェッチするPDB IDを手動で指定します。
    rna_pdb_ids = load_pdb_ids()

    # 指定されたディレクトリにPDBファイルをダウンロード
    for pdb_id in rna_pdb_ids:
        pdbl.retrieve_pdb_file(pdb_id, pdir='/work/gs58/s58007/rna_GWclustering/data/pdb', file_format='pdb')

# 関数を実行してRNAのPDBファイルをダウンロード
fetch_all_rna_structures()