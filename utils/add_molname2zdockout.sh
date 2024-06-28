# export input="R0254/zdock/pdb/zdock.S_000239-top1.pdb"
export input="R0254/pymol_converted/sample239.pdb"
export output="R0254/zdock.moladded.000239.pdb"

# cut で 10, 12, 13 列目を消す。そして末尾に element を追加
# 1列目が ATOM のところだけ
cat $input | \
awk '{
    if ($1 != "ATOM") {
        print $0;
        next;
    }
    atom_name = $3;  # 原子名が3番目の列にあると仮定
    element = substr(atom_name, 1, 1);  # 原子名の最初の1文字を抽出

    # 元素記号が2文字の場合の処理（例：Cl, Brなど）
    if (atom_name ~ /^[A-Z][a-z]/) {
        element = substr(atom_name, 1, 2);
    }
    # タブで区切って element を末尾に追加
    print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9  "\t"  element "\t" $11;
}' > $output