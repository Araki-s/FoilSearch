from tqdm import tqdm
import numpy as np
import warnings
from scipy.optimize import curve_fit
import os
import math
from time import sleep

reference = 'reference/'
output = 'OUTPUT/'
warnings.simplefilter('ignore')

#便利機能###########################################################################################################
def Read_tab_splited(x):
    """
    タブ区切りのファイルからデータを読み取る関数
    引数：ファイルのパス
    返り値：行ごとのリストのリスト
    """
    tab_splited = []
    f = open(x)
    line = f.readline()
    while line:
        tab_splited.append(line.rstrip().split("\t"))
        line = f.readline()
    f.close
    return tab_splited 

def Read_n_splited(x):
    """
    改行区切りのファイルからデータを読み取る関数
    引数：ファイルのパス
    返り値：行ごとのリストのリスト
    """
    n_splited = []
    f = open(x)
    data1 = f.read()
    f.close()
    lines1 = data1.split("\n") 
    for line in lines1:
        n_splited.append(line+'  ')
    return n_splited

def numeric(a):
    """
    JENDLから読み取った断面積値、エネルギー値の'1.000+6'を'1.000E+6'に変える関数
    引数：JENDLから読み取った文字列(ex:1.000+6)
    返り値:'1.000E+6'
    """
    k = ""
    if a[len(a)-2] == "+" or a[len(a)-2] == "-":
        for i in range(len(a)):
            if i == len(a)-2:
                k = k + "E" + a[i]
            else:
                k = k + a[i]
        return k
    if a[len(a)-3] == "+" or a[len(a)-3] == "-":
        for i in range(len(a)):
            if i == len(a)-3:
                k = k + "E" + a[i]
            else:
                k = k + a[i]
        return k
    else:
        return a
    
def Calculate_interpolation(x,y,z):
    """
    データの内挿を計算する関数
    x,y:内挿をするデータのリスト(数値)
    ☆ len(x)=len(y)出ないと動かないです。
    z:内挿するｘ座標（数値）
    """
    if len(x) != len(y):
        print('適切に引数を設定してください。FOIL.pyを参照のこと')
        return 'error in Calculate_interpolation'
    result = []
    for i in range(len(z)):
        done = False
        if z[i] < x[0]:
            done = True
            interpolation = ((y[1] - y[0])/(x[1] - x[0]))*z[i] + ((y[0]*x[1] - y[1]*x[0])/(x[1] - x[0]))
        elif z[i] > x[len(x)-1]:
            done = True
            interpolation = ((y[len(x)-1] - y[len(x)-2])/(x[len(x)-1] - x[len(x)-2]))*z[i] + ((y[len(x)-2]*x[len(x)-1] - y[len(x)-1]*x[len(x)-2])/(x[len(x)-1] - x[len(x)-2]))
        else:
            for j in range(len(x)-1):
                if z[i] >= x[j] and z[i] <= x[j-1]:
                    done = True
                    interpolation = ((y[j] - y[j-1])/(x[j] - x[j-1]))*z[i] + ((y[j-1]*x[j] - y[j]*x[j-1])/(x[j] - x[j-1]))
            if not done:
                 print('適切に引数を設定してください。FOIL.pyを参照のこと')
                 return 'error in Calculate_interpolation'
        result.append(interpolation)        
    return result

def Read_comma_splited(x):
    """
    カンマ区切りのファイルからデータを読み取る関数
    引数：ファイルのパス
    返り値：行ごとのリストのリスト
    """
    comma_splited = []
    f = open(x)
    line = f.readline()
    while line:
        comma_splited.append(line.rstrip().split(','))
        line = f.readline()
    f.close
    return comma_splited

def write_file_nested_list(path,write_list):
    with open(path, mode = "w") as f:
        for i,name_i in enumerate(write_list):
            for ii, name_ii in enumerate(name_i):
                if ii != len(name_i)-1:
                    f.write(str(name_ii) + '\t')
                else:
                    f.write(str(name_ii))
            if i != len(write_list)-1:
                f.write('\n')    

#Phase1#############################################################################################################
def Phase_1(efficiency_file):
    """
    箔の最適化の第一段階をまとめた関数
    標準線源での実験結果をもとに、Ge検出器の効率関数f(E)を求める
    """
    print('Phase_1開始')
    Curve_fit(Read_tab_splited(reference + efficiency_file))
    print('Phase_1終了')

def Curve_fit(y):
    """
    検出効率を最小二乗法で求める関数
    パラメータaとEをCurve_fit.txtにタブ区切りで出力
    引数：[エネルギー、効率]のリスト
    """
    
    xdata = []
    ydata = []
    for i in y:
        xdata.append(float(i[0]))
        ydata.append(float(i[1]))
    x = []
    for i in range(50,1300):
        x.append(i)
    ##### 最適化するy=f(x)の係数の初期値
    prameter_initial = np.array([1, 1, 1, 1, 1, 1, 1])

    ##### 最適化の計算[2]
    popt, pcov = curve_fit(func, xdata, ydata, p0= prameter_initial)

    ##### 最適化後のy=f(x)の関数
    y = func(x, *popt)
    
    '''
    ##### 生データと最適化後の関数をグラフにプロット
    # 元の生データを点プロット
    plt.scatter(xdata, ydata, c='blue', label='raw data')
    # 最適化後のフィット関数を線でプロット
    plt.plot(x, y, 'r-',label = "fitting")
    
    ##### グラフ表示のオプション設定 #####
    plt.xlabel('Neutron Energy [KeV]', fontsize=18)     # x軸の名称とフォントサイズ
    plt.ylabel('Efficiency', fontsize=18)               # y軸の名称とフォントサイズ
    plt.yscale("log")                                   # y軸のスケールをログに指定
    plt.xscale("log")  
    plt.legend(loc='upper right')                       # ラベルを右上に記載
    plt.tight_layout()                                  # タイトルを重ねない
    #plt.show()                                         # グラフをプロット
    plt.savefig( output + 'Curve_fit.png')              # 画像をファイルで保存
    '''
    with open(output + 'Curve_fit.txt' , mode = "w") as f:
        for i, name in enumerate(popt):
            if i != len(popt)-1:
                f.write(str(name)+',')
            else:
                f.write(str(name))

def func(x, a0, a1, a2, a3 , a4, a5, E):
    """
    検出効率をフィッティングするための関数
    """
    return a0 + a1*np.log(x/E) + a2*(np.log(x/E))**2 + a3*(np.log(x/E))**3 + a4*(np.log(x/E))**4 + a5*(np.log(x/E))**5

#Phase2##############################################################################################################
def Phase_2(NIST_file,JENDL,density_file,Gammas_file,attenuation_folder,Threshold_MIN,Threshold_MAX):
    print('Phase_2開始')
    sleep(1)
    if os.path.exists(output + 'result_list_' + str(Threshold_MIN) + '_' + str(Threshold_MAX) + '.txt' ):
        pass
    else:
        Make_result_list(NIST_file,JENDL,density_file,Gammas_file,attenuation_folder,Threshold_MIN,Threshold_MAX)
    print('Phase_2終了')
    sleep(1)

def Read_NIST(NIST_file):
    """
    NIST.txtからデータを読み取る関数
    NIST.txt→核種ごとの情報(原子番号、元素記号、質量数,etc)ファイル
    """
    
    natural_element_nist_list = []#自然界に存在する核種の元素記号
    natural_comp_list = []#自然界に存在する核種の存在比
    natural_ramass_list = []#自然界に存在する核種の相対原子質量
    natural_massnum_list = []#自然界に存在する核種の原子番号
    natural_saw_list = []#自然界に存在する核種の平均質量数
    element_nist_list = []#すべての元素記号
    massnum_list = []#すべての原子番号
    nist_list = Read_n_splited(reference + NIST_file)

    for i in range(math.floor(len(nist_list)/8)):
        if len(nist_list[8*i+4]) != 25:#Isotopic Compositionのデータがあるとき（存在比が有意）
            if nist_list[8*i+2][15] == " ":#質量数が１桁の時
                natural_element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17]+"-"+nist_list[8*i+2][14].rjust(3))
            elif nist_list[8*i+2][16] == " ":#質量数が2桁の時
                natural_element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17]+"-"+(nist_list[8*i+2][14]+nist_list[8*i+2][15]).rjust(3))
            else:#質量数が３桁の時
                natural_element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17]+"-"+nist_list[8*i+2][14]+nist_list[8*i+2][15]+nist_list[8*i+2][16])

            #同位体存在比を抽出--------------------
            s = nist_list[8*i+4].split()
            comp = ""
            for k in range(len(s[3])):
                if s[3][k] == "(":
                    break
                comp = comp + s[3][k]
            natural_comp_list.append(comp.ljust(10))

            #相対原子質量を抽出--------------------
            s = nist_list[8*i+3].split()
            ramass = ""
            for l in range(len(s[4])):
                if s[4][l] == "(":
                    break
                ramass = ramass + s[4][l]
            natural_ramass_list.append(ramass)

            #原子番号を抽出---------------------
            natural_massnum_list.append(int(nist_list[8*i][16]+nist_list[8*i][17]+nist_list[8*i][18]))

            #平均原子量を抽出
            s = nist_list[8*i+5].split()
            saw = ""
            for j in range(len(s[4])):
                if s[4][j] == "(":
                    break
                if s[4][j] != "[" and s[4][j] != "]":
                    saw = saw + s[4][j]
            if len(saw.split(",")) == 2:
                d = saw.split(",")
                saw = (float(d[0]) + float(d[1]))/2
            natural_saw_list.append(saw)


    #すべての元素を抽出------------------------------------------------------------------------------
    for i in range(math.floor(len(nist_list)/8)):
        element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17])

        #原子番号を抽出---------------------
        massnum_list.append(int(nist_list[8*i][16]+nist_list[8*i][17]+nist_list[8*i][18]))    

def Make_result_list(NIST_file,JENDL,density_file,Gammas_file,attenuation_folder,Threshold_MIN,Threshold_MAX):
    #NIST.txtからの読み取り部分___________________________________________________________________
    reaction_list = []
    result_list = []
    natural_element_nist_list = []#自然界に存在する核種の元素記号
    natural_comp_list = []#自然界に存在する核種の存在比
    natural_ramass_list = []#自然界に存在する核種の相対原子質量
    natural_massnum_list = []#自然界に存在する核種の原子番号
    natural_saw_list = []#自然界に存在する核種の平均質量数
    element_nist_list = []#すべての元素記号
    massnum_list = []#すべての原子番号
    nist_list = Read_n_splited(reference + NIST_file)

    for i in range(math.floor(len(nist_list)/8)):
        if len(nist_list[8*i+4]) != 25:#Isotopic Compositionのデータがあるとき（存在比が有意）
            if nist_list[8*i+2][15] == " ":#質量数が１桁の時
                natural_element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17]+"-"+nist_list[8*i+2][14].rjust(3))
            elif nist_list[8*i+2][16] == " ":#質量数が2桁の時
                natural_element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17]+"-"+(nist_list[8*i+2][14]+nist_list[8*i+2][15]).rjust(3))
            else:#質量数が３桁の時
                natural_element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17]+"-"+nist_list[8*i+2][14]+nist_list[8*i+2][15]+nist_list[8*i+2][16])

            #同位体存在比を抽出--------------------
            s = nist_list[8*i+4].split()
            comp = ""
            for k in range(len(s[3])):
                if s[3][k] == "(":
                    break
                comp = comp + s[3][k]
            natural_comp_list.append(comp.ljust(10))

            #相対原子質量を抽出--------------------
            s = nist_list[8*i+3].split()
            ramass = ""
            for l in range(len(s[4])):
                if s[4][l] == "(":
                    break
                ramass = ramass + s[4][l]
            natural_ramass_list.append(ramass)

            #原子番号を抽出---------------------
            natural_massnum_list.append(int(nist_list[8*i][16]+nist_list[8*i][17]+nist_list[8*i][18]))

            #平均原子量を抽出
            s = nist_list[8*i+5].split()
            saw = ""
            for j in range(len(s[4])):
                if s[4][j] == "(":
                    break
                if s[4][j] != "[" and s[4][j] != "]":
                    saw = saw + s[4][j]
            if len(saw.split(",")) == 2:
                d = saw.split(",")
                saw = (float(d[0]) + float(d[1]))/2
            natural_saw_list.append(saw)

    #すべての元素を抽出----------------------------
    for i in range(math.floor(len(nist_list)/8)):
        element_nist_list.append(nist_list[8*i+1][16]+nist_list[8*i+1][17])

        #原子番号を抽出---------------------
        massnum_list.append(int(nist_list[8*i][16]+nist_list[8*i][17]+nist_list[8*i][18]))    

    #JENDLからの読み取り部分________________________________________________________________________
    cross_section_x_list = []#ｘとｙに分けた断面積スペクトル
    cross_section_y_list = []#ｘとｙに分けた断面積スペクトル
                                
    #ファイル名を取得し、ファイルごとに情報を抜き出す
    files = os.listdir(reference + JENDL)
    for t in tqdm(files):
        text_list = []#１行ずつテキストを格納する配列
        MF10_existence = 0#MF=10があるかどうかを識別する定数（0:なし　１:あり）
        MT_list = []#MF=10のMTが格納される配列 MF=10 MT=AAA (X,X)
        MF10Q0 = 0#基底順位のQ値を代入する変数
        MF10Q = 0#Q値を代入する変数
        MF10Qline_list = []#放射化後の核種ごとのQ値が書かれている行が格納される配列
        reaction_per_element_list = []#１つの核種のすべての反応をまとめたリスト[核種、反応、順位、閾値、娘核の原子番号、娘核の質量数、娘核の核種、存在比、平均質量]のリスト
        cross_section_per_element_list = []#反応ごとの断面積スペクトル
        ll=["(n,He3)","(n,t)","(n,2na)","(n,2p)","(n,d)","(n,2n)","(n,np)","(n,p)","(n,pa)","(n,na)","(n,n')","(n,a)","(n,nd)","(n,2np)","(n,3n)","(n,nt)"]#自分が書いた反応
        num_lines = len(open(reference + JENDL + t).readlines())#最終行を取得
        lines = Read_n_splited(reference + JENDL + t)
        
        n=0
        for line in lines:
            n += 1
            text_list.append(line)
            
            #MF=10の記載があるかどうかを識別
            if line.startswith("MF=10"):
                MF10_existence = 1
                
            #MF=10がある時
            if MF10_existence == 1:
                if line.startswith("  MT"):#行が MTから始まるとき、その行の必要な情報をMT_listに追加
                    MT_list.append("MF=10 "+ line[2]+line[3]+line[4]+line[5]+line[6]+line[7]+line[8]+line[9]+line[10]+line[11]+line[12]+line[13]+line[14]+line[15]+line[16]\
                                   +line[17]+line[18]+line[19]+line[20]+line[21]+line[22]+line[23]+line[24]+line[25])
            if n != num_lines + 1:
                for i in range(len(MT_list)):
                    if line[70]+line[71]+line[72]+line[73]+line[74] == "10" + MT_list[i][9]+MT_list[i][10]+MT_list[i][11]:#71～74列目が10MTに一致する行（MF=10,MT=〇）
                        if line[78] + line[79] == " 2":#行末が２に一致する行（基底エネルギーが書いてある行）
                            MF10Q0 = line[0]+line[1]+line[2]+line[3]+line[4]+line[5]+line[6]+line[7]+line[8]+line[9]+line[10]
                        if line[0]+line[1]+line[2]+line[3]+line[4]+line[5]+line[6]+line[7]+line[8]+line[9]+line[10] == MF10Q0:#Q値が書いてある行
                            MF10Qline_list.append(n-1)
                            
        #１つの核種の反応ごとのQ値などを読み取る
        #１つの核種の反応ごとにforを回している----------------------------------------------------------------
        for i in range(len(MF10Qline_list)):
            data_list = []#forループで回っている反応ごとの情報が追加されるリスト
            cross_data_list = []#forループで回っている反応ごとの断面積スペクトルが追加されるリスト
            k = MF10Qline_list[i]#kにQ値が書いてある行を代入
            
            data_list.append(text_list[0][16]+text_list[0][17]+text_list[0][18]+text_list[0][19]+text_list[0][20]+text_list[0][21])#放射化前の核種
            for j in range(len(MT_list)):
                if MT_list[j][9]+MT_list[j][10]+MT_list[j][11] == text_list[k][72]+text_list[k][73]+text_list[k][74]:
                    data_list.append(MT_list[j])#反応の種類
                    
            #Q値
            MF10Q = numeric(text_list[k][11]+text_list[k][12]+text_list[k][13]+text_list[k][14]+text_list[k][15]+text_list[k][16]\
                            +text_list[k][17]+text_list[k][18]+text_list[k][19]+text_list[k][20]+text_list[k][21])
            
            #グランドステイトのQ値
            MF10Q0 = numeric(text_list[k][0]+text_list[k][1]+text_list[k][2]+text_list[k][3]+text_list[k][4]+text_list[k][5]\
                             +text_list[k][6]+text_list[k][7]+text_list[k][8]+text_list[k][9]+text_list[k][10])
            
            #グランドステイトからの差分
            MF10dQ = (float(MF10Q0) - float(MF10Q))/1000
            data_list.append(MF10dQ)
            
            #Q値を追加-------------------------------------------------
            if float(MF10Q) < 0:#Q値が負の時（＝閾値がある）
                #閾値
                MF10Threshold = numeric(text_list[k+2][1]+text_list[k+2][2]+text_list[k+2][3]+text_list[k+2][4]+text_list[k+2][5]+text_list[k+2][6]\
                                        +text_list[k+2][7]+text_list[k+2][8]+text_list[k+2][9]+text_list[k+2][10]+text_list[k+2][11])
                data_list.append(MF10Threshold)
            else:#Q値が正の時（＝閾値がなく、14MeV中性子に対する断面積の情報はいらない）
                data_list.append("0")#閾値がないので閾値の欄に0を追加
                
            #断面積スペクトルを取得-------------------------------------------------
            if i != len(MF10Qline_list)-1:#最終行まで
                #反応ごとのforが対象にしている反応を１列ずつ読む----
                for u in range(MF10Qline_list[i]+2,MF10Qline_list[i+1]):
                    if text_list[u][76]+text_list[u][77]+text_list[u][78]+text_list[u][79] != "9999":
                        if text_list[u][78]+text_list[u][79] != " 1":
                            for a in range(len(text_list[u].split())):
                                if len(text_list[u].split()[a]) >= 8 and ("+" in text_list[u].split()[a]  or "-" in text_list[u].split()[a] ):
                                    cross_data_list.append(numeric(text_list[u].split()[a][0]+text_list[u].split()[a][1]+text_list[u].split()[a][2]\
                                                                   +text_list[u].split()[a][3]+text_list[u].split()[a][4]+text_list[u].split()[a][5]\
                                                                   +text_list[u].split()[a][6]+text_list[u].split()[a][7]+text_list[u].split()[a][8]+text_list[u].split()[a][9]))
            else:
                for u in range(MF10Qline_list[i]+2,num_lines):
                    if text_list[u][76]+text_list[u][77]+text_list[u][78]+text_list[u][79] == "9999":
                        break
                    if text_list[u][76]+text_list[u][77]+text_list[u][78]+text_list[u][79] != "9999":
                        if text_list[u][78]+text_list[u][79] != " 1":
                            for a in range(len(text_list[u].split())):
                                if len(text_list[u].split()[a]) >= 8 and ("+" in text_list[u].split()[a]  or "-" in text_list[u].split()[a] ):
                                    cross_data_list.append(numeric(text_list[u].split()[a][0]+text_list[u].split()[a][1]+text_list[u].split()[a][2]+text_list[u].split()[a][3]\
                                                                   +text_list[u].split()[a][4]+text_list[u].split()[a][5]+text_list[u].split()[a][6]+text_list[u].split()[a][7]\
                                                                   +text_list[u].split()[a][8]+text_list[u].split()[a][9]))  
            cross_section_per_element_list.append(cross_data_list)#反応ごとの断面積スペクトルを核種でまとめる
            reaction_per_element_list.append(data_list)#それぞれの反応の情報を核種でまとめる
        
        if len(reaction_per_element_list) != 0:
            for i in range(len(reaction_per_element_list)):
                for ii in range(len(natural_element_nist_list)):
                    if natural_element_nist_list[ii] == reaction_per_element_list[i][0]:#自然界に存在するとき
                        
                        #閾値の条件を満たすとき
                        if float(reaction_per_element_list[i][3]) <= Threshold_MAX and float(reaction_per_element_list[i][3]) >= Threshold_MIN:
                            
                            ##じぶんで書いた以下の反応のif分に漏れがないか調べるためのprint----------------------------------------
                            ext = 0
                            for iii in ll:
                                if iii in reaction_per_element_list[i][1]:
                                    ext =1
                            if ext != 1:
                                print(reaction_per_element_list[i],"この反応が漏れています")#じぶんで書いた以下の反応のif分に漏れがないか調べるためのprint  
                                
                            if  "(n,He3)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-2:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-2)#原子番号が2減る
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)#質量数が2減る
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                    
                            if  "(n,t)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-1:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-1)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                            
                            if  "(n,2na)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-2:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-2)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-5)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                            
                            if  "(n,2p)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-2:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-2)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,d)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-1:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-1)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,2n)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii])
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,np)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-1:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-1)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,p)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-1:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-1)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5]))
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                            
                            if  "(n,pa)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-3:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-3)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-4)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,na)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-2:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-2)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-4)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,n')" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii])
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5]))
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,a)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-2:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-2)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-3)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,nd)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-1:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-1)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,2np)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-1:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-1)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
    
                            if  "(n,3n)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii])
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
                                        
                            if  "(n,nt)" in reaction_per_element_list[i][1]:
                                for j in range(len(massnum_list)):
                                    if massnum_list[j] == natural_massnum_list[ii]-1:
                                        reaction_per_element_list[i].append(natural_massnum_list[ii]-1)
                                        reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-3)
                                        reaction_per_element_list[i].append(element_nist_list[j].rstrip())
                                        reaction_per_element_list[i].append(natural_comp_list[ii].rstrip())
                                        reaction_per_element_list[i].append(natural_saw_list[ii])
                                        break
    
                            #reaction_per_element_listが条件を満たすときその断面積スペクトルを取得-------------------
                            x_list = []
                            y_list = []
                            
                            for iii in range(int(len(cross_section_per_element_list[i])/2)):
                                x_list.append(float(cross_section_per_element_list[i][2*iii]))
                                y_list.append(float(cross_section_per_element_list[i][2*iii+1]))
                
                            #断面積スペクトルを保存-----------------------------------
                            cross_section_x_list.append(x_list)
                            cross_section_y_list.append(y_list)
                            reaction_list.append(reaction_per_element_list[i])
    density_list = Read_comma_splited(reference + density_file)
    Gammas_list = Read_comma_splited(reference + Gammas_file)
    cross_section_per_gammas_x_list = []
    cross_section_per_gammas_y_list = []
    
    for i in tqdm(range(len(reaction_list))):
        errors = []
        k = 0 #jendlから読み取った放射化後の順位のデータが自分で作ったgammmas.csvにあるかどうかを判別するための変数     
        #reaction_listに密度を追加-------------------------------------------------------------------------------------------------
        pp = 0#密度データがdensity.csvにあるかどうかを判別するための変数
        for ii in range(len(density_list)):
            if density_list[ii][0].ljust(2) == reaction_list[i][0][0]+reaction_list[i][0][1]:
                reaction_list[i].append(density_list[ii][1])
                pp = 1
        if pp == 0:
            errors.append("density.csvに",reaction_list[i][0][0]+reaction_list[i][0][1],"を追加してください")
            reaction_list[i].append("密度データがありません")
            
        #ガンマ線のエネルギー、放射率、半減期を読み取る--------------------------------------------------------------------------------
        for ii in range(len(Gammas_list)):
            
            if "x" not in Gammas_list[ii][2] and "y" not in Gammas_list[ii][2] and "x" not in Gammas_list[ii][3]:
                
                if float(reaction_list[i][4]) ==float(Gammas_list[ii][0]) and float(reaction_list[i][5]) == float(Gammas_list[ii][1]) \
                and (round(float(reaction_list[i][2]),5) < float(Gammas_list[ii][2])+3 and round(float(reaction_list[i][2]),5) > float(Gammas_list[ii][2])-3):#reaction_list(reaction_list)がgammmas.csvと一致するとき    
                    a = []
                    k=1    
                    #reaction_listの内容をaリストに追加
                    for t in range(len(reaction_list[i])):
                        a.append(reaction_list[i][t])
    
                    if float(Gammas_list[ii][3]) != 0:#ガンマ線のエネルギー!=0のとき＝ステイブルではないとき
                        a.append(Gammas_list[ii][3])#ガンマ線のエネルギー
                        a.append(Gammas_list[ii][4])#放出率
                        a.append(Gammas_list[ii][5])#半減期
    
    
                        exsistence_attenuation_file = 0#吸収係数のファイルがあるかどうかをみる変数
                        Attenuation_files = os.listdir(reference + attenuation_folder)
                        for q in Attenuation_files:
                            if q == (reaction_list[i][0][0]+reaction_list[i][0][1]).rstrip()+".txt":
                                exsistence_attenuation_file = 1
                        if exsistence_attenuation_file == 0:
                            print("Attenuationに",(reaction_list[i][0][0]+reaction_list[i][0][1]).rstrip(),"がありません。追加してください。")
                            a.append("Attenuationのデータがありません")
                        else:        
                            #質量吸収係数を読み取る、a2に入れる-----------------------------------------------------------------------
                            a2 = Read_tab_splited(reference + attenuation_folder +'/%s.txt'%(reaction_list[i][0][0]+reaction_list[i][0][1]).rstrip())
                            for iii in range(len(a2[0])):
                                a2[0][iii] = float(a2[0][iii])
                                a2[1][iii] = float(a2[1][iii])
                            #質量吸収係数を算出するためにガンマ線のエネルギーが質量吸収係数のデータのどの座標にあるかを調べる--------------
                            attenuation = Calculate_interpolation(a2[0], a2[1], [float(Gammas_list[ii][3])*10**(-3)])[0]
                        a.append(attenuation)

                        efficiency = func(float(Gammas_list[ii][3]),*[float(x) for x in Read_comma_splited(output+'Curve_fit.txt')[0]])
                        a.append(efficiency)
    
                    elif float(Gammas_list[ii][3]) == 0:
                        a.append(Gammas_list[ii][3])
                        a.append(Gammas_list[ii][4])
                        a.append(Gammas_list[ii][5])
                        a.append(0)
                        a.append(0)
    
                    result_list.append(a)
                    cross_section_per_gammas_x_list.append(cross_section_x_list[i])
                    cross_section_per_gammas_y_list.append(cross_section_y_list[i])
    
        #if k == 0:#jendlから読み取ったデータがgammmas.csvにない
            #print(reaction_list[i][4],reaction_list[i][5],reaction_list[i][6],reaction_list[i][2],reaction_list[i][1])

    write_file_nested_list(output + 'result_list_' + str(Threshold_MIN) + '_' + str(Threshold_MAX) + '.txt',result_list)
    write_file_nested_list(output + 'cross_section_x' + str(Threshold_MIN) + '_' + str(Threshold_MAX) + '.txt',cross_section_per_gammas_x_list)
    write_file_nested_list(output + 'cross_section_y' + str(Threshold_MIN) + '_' + str(Threshold_MAX) + '.txt',cross_section_per_gammas_y_list)

def read_results(Threshold_MIN,Threshold_MAX):
    file_list =[output + 'result_list_' + str(Threshold_MIN) + '_' + str(Threshold_MAX) + '.txt',output + 'cross_section_x' + str(Threshold_MIN) + '_' + str(Threshold_MAX) + '.txt',output + 'cross_section_y' + str(Threshold_MIN) + '_' + str(Threshold_MAX) + '.txt']
    result_list = []
    cross_section_x = []
    cross_section_y = []
    for i ,file in enumerate(file_list):
        f = open(file)
        data = f.readlines()
        if i == 0:
            for j in data:
                result_list.append(j.split('\t'))
        elif i == 1:
            for j in data:
                cross_section_x.append(j.split('\t'))
        elif i == 2:
            for j in data:
                cross_section_y.append(j.split('\t'))
        f.close()
    for i in range(len(result_list)):
        result_list[i][3] = float(result_list[i][3])#閾値
        result_list[i][7] = float(result_list[i][7])#存在比
        result_list[i][8] = float(result_list[i][8])#平均原子量
        result_list[i][9] = float(result_list[i][9])#密度
        result_list[i][10] = float(result_list[i][10])#ガンマ線のエネルギー
        result_list[i][13] = float(result_list[i][13])#質量吸収係数
        result_list[i][14] = float(result_list[i][14])#検出効率
        #↓放出比
        if result_list[i][11] == '*':
            result_list[i][11] = 0.0
        else:
            result_list[i][11] = float(result_list[i][11].replace('<','').replace('>','').replace('~','').replace('*',''))/100
        #↓半減期↓
        temp = ''
        if len(result_list[i][12]) > 1:
            for ii in range(len(result_list[i][12])-1):
                temp = temp + result_list[i][12][ii]
        else:
            temp = result_list[i][12]
        temp = float(temp)
        if result_list[i][12][len(result_list[i][12])-1] == "s":
            temp = temp
        elif result_list[i][12][len(result_list[i][12])-1] == "m":
            temp = temp*60
        elif result_list[i][12][len(result_list[i][12])-1] == "h":
            temp = temp*60*60
        elif result_list[i][12][len(result_list[i][12])-1] == "d":
            temp = temp*60*60*24
        elif result_list[i][12][len(result_list[i][12])-1] == "y":
            temp = temp*60*60*24*365
        result_list[i][12] = float(temp)
    for i in range(len(cross_section_x)):
        for j in range(len(cross_section_x[i])):
            cross_section_x[i][j] = float(cross_section_x[i][j])
            cross_section_y[i][j] = float(cross_section_y[i][j])
        
    return result_list, cross_section_x, cross_section_y

