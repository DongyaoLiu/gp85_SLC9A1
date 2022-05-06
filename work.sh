python3 /Users/liudongyao/Downloads/NWAFU/ALV/parse_GenBank_gb.py ./sequence.gb.txt faa >01.all.alv.faa
python3 /Users/liudongyao/Downloads/NWAFU/ALV/parse_GenBank_gb.py ./sequence.gb.txt fa >01.all.alv.fa
python3 /Users/liudongyao/Downloads/NWAFU/ALV/parse_GenBank_gb.py ./sequence.gb.txt annotation | cut -d">" -f2 >01.all.alv.annotation.txt
