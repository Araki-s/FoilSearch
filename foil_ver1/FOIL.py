from tqdm import tqdm as tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import math

class Optimum_foil:
    def __init__(self,input_file,output_file,efficiency_file,NIST_file,JENDL_file,density_file,Gammas_file,attenuation_file,spectrum_file,\
                 Threshold_MIN,Threshold_MAX,Threshold_MIN_calc,Threshold_MAX_calc,Threshold_MIN_plot,Threshold_MAX_plot,day,tc,\
                     plot_activation_cross_section,plot_lower_limit):

        
        self.input_file = input_file
        self.output_file = output_file
        self.efficiency_file = efficiency_file
        self.NIST_file = NIST_file
        self.JENDL_file = JENDL_file
        self.density_file = density_file
        self.Gammas_file = Gammas_file
        self.attenuation_file = attenuation_file
        self.spectrum_file = spectrum_file
        self.Threshold_MIN = Threshold_MIN
        self.Threshold_MAX = Threshold_MAX
        self.Threshold_MIN_calc = Threshold_MIN_calc
        self.Threshold_MAX_calc = Threshold_MAX_calc
        self.Threshold_MIN_plot = Threshold_MIN_plot
        self.Threshold_MAX_plot = Threshold_MAX_plot
        self.day = day
        self.tc = tc
        self.plot_activation_cross_section = plot_activation_cross_section
        self.plot_lower_limit = plot_lower_limit
        

        self.reaction_list = []
        self.spectrum_x_list = []
        self.spectrum_y_list = []
        self.cross_section_Nb_list = []
        self.result_list = []
        self.sigma_reaction_Nb = 0
        

        
    def Phase_1(self):
        """
        箔の最適化の第一段階をまとめた関数
        """
        self.Curve_fit(self.Read_tab_splited(self.input_file + self.efficiency_file))
             
    def Phase_2(self):
        if os.path.exists(self.output_file + '/Phase2/reaction_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt' ):
            print('Phase2で作成するファイルがすでに在します')
        else:
            self.Read_NIST()
            self.Make_reaction_list()
            
    def Phase_3(self):
        if os.path.exists(self.output_file + '/Phase3/result_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt' ):
            print('Phase3で作成するファイルがすでに在します')
        else:
            self.Make_result_list()
    
    def write_file_nested_list(self,path,write_list):
        with open(path, mode = "w") as f:
            for i,name_i in enumerate(write_list):
                for ii, name_ii in enumerate(name_i):
                    if ii != len(name_i)-1:
                        f.write(str(name_ii) + '\t')
                    else:
                        f.write(str(name_ii))
                if i != len(write_list)-1:
                    f.write('\n')


    def write_file_non_nested_list(self,path,write_list):
        with open(path, mode = "w") as f:
            for k,name in enumerate(write_list):
                if k != len(write_list)-1:
                    f.write(str(name))
                    f.write('\n')
                else:
                    f.write(str(name))
                    
    def numeric(self,a):
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
    
    def get_unique_list(self,seq):
        seen = []
        return [x for x in seq if x not in seen and not seen.append(x)]
        
    def Read_tab_splited(self, x):
        """
        タブ区切りのファイルからデータを読み取る関数
        引数：ファイルのパス
        返り値：行ごとのリストのリスト
        """
        self.tab_splited = []
        f = open(x)
        line = f.readline()
        while line:
            self.tab_splited.append(line.rstrip().split("\t"))
            line = f.readline()
        f.close
        return self.tab_splited 
    
    def Read_n_splited(self ,x):
        """
        改行区切りのファイルからデータを読み取る関数
        引数：ファイルのパス
        返り値：行ごとのリストのリスト
        """
        self.n_splited = []
        f = open(x)
        data1 = f.read()
        f.close()
        lines1 = data1.split("\n") 
        for line in lines1:
            self.n_splited.append(line+'  ')
        return self.n_splited
    
    def Read_comma_splited(self,x):
        """
        カンマ区切りのファイルからデータを読み取る関数
        引数：ファイルのパス
        返り値：行ごとのリストのリスト
        """
        self.comma_splited = []
        f = open(x)
        line = f.readline()
        while line:
            self.comma_splited.append(line.rstrip().split(','))
            line = f.readline()
        f.close
        return self.comma_splited

    def func(self,x, a0, a1, a2, a3 , a4, a5, E):
        """
        検出効率をフィッティングするための関数
        """
        return a0 + a1*np.log(x/E) + a2*(np.log(x/E))**2 + a3*(np.log(x/E))**3 + a4*(np.log(x/E))**4 + a5*(np.log(x/E))**5
        
    def Curve_fit(self, y):
        """
        検出効率を最小二乗法で求める関数
        パラメータaとEをCurve_fit.txtにタブ区切りで出力
        引数：[エネルギー、効率]のリスト
        """
        self.xdata = []
        self.ydata = []
        for i in y:
            self.xdata.append(float(i[0]))
            self.ydata.append(float(i[1]))
        x = []
        for i in range(50,1300):
            x.append(i)


        ##### 最適化するy=f(x)の係数の初期値
        prameter_initial = np.array([1, 1, 1, 1, 1, 1, 1])

        ##### 最適化の計算[2]
        popt, pcov = curve_fit(self.func, self.xdata, self.ydata, p0= prameter_initial)
        print ("parameter ->", popt)

        ##### 最適化後のy=f(x)の関数
        y = self.func(x, *popt)

        ##### 生データと最適化後の関数をグラフにプロット
        # 元の生データを点プロット
        plt.scatter(self.xdata, self.ydata, c='blue', label='raw data')
        # 最適化後のフィット関数を線でプロット
        plt.plot(x, y, 'r-',label = "fitting")

        ##### グラフ表示のオプション設定 #####
        plt.xlabel('Neutron Energy [KeV]', fontsize=18)     # x軸の名称とフォントサイズ
        plt.ylabel('Efficiency', fontsize=18)     # y軸の名称とフォントサイズ
        plt.yscale("log")                # y軸のスケールをログに指定
        plt.xscale("log")  
        plt.legend(loc='upper right')    # ラベルを右上に記載
        plt.tight_layout()   # タイトルを重ねない
        #plt.show()                       # グラフをプロット
        plt.savefig(self.output_file + '/Curve_fit.png')  # 画像をファイルで保存
        
        with open(self.output_file + '/Curve_fit.txt' , mode = "w") as f:
            for i, name in enumerate(popt):
                if i != len(popt)-1:
                    f.write(str(name)+',')
                else:
                    f.write(str(name))
    
    def Calculate_interpolation(self,x,y,z,mode):
        """
        データの内挿を計算する関数
        x,y:内挿をするデータのリスト(数値)
        z:内挿するｘ座標（数値）
        mode:1=外挿していることを通知、2=外挿していることを通知しない
        """
        interpolation = 0
        interpolation_line = 0
        if mode == 1:
            if float(z) > float(x[0]):
                print('下限で外挿しています')
            elif float(z) > float(x[len(x)-1]):
                print('上限で外挿しています')
                
        for i in range(len(x)):
            if float(z) > float(x[i]):
                interpolation_line = i
                
        
        if interpolation_line == len(x)-1:
            x1 = float(x[interpolation_line-1])
            y1 = float(y[interpolation_line-1])
            x2 = float(x[interpolation_line])
            y2 = float(y[interpolation_line])
            x = float(z)

        else:
            x1 = float(x[interpolation_line])
            y1 = float(y[interpolation_line])
            x2 = float(x[interpolation_line+1])
            y2 = float(y[interpolation_line+1])
            x = float(z)
            
        interpolation = (y2-y1)/(x2-x1)*(x-x1)+y1
        return interpolation
       
    def plot_isomer(self):
        result_list = self.Read_tab_splited(self.output_file + '/Phase3/result_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_x_list = self.Read_tab_splited(self.output_file + '/Phase3/cross_section_x' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_y_list = self.Read_tab_splited(self.output_file + '/Phase3/cross_section_y' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        isomer_list = []
        
        for i in tqdm(range(len(result_list))):
            if float(result_list[i][2]) != 0 and max([float(k) for k in cross_section_y_list[i]]) > 1:
                plt.plot([float(k) for k in cross_section_x_list[i]],[float(k) for k in cross_section_y_list[i]])
                isomer_list.append(result_list[i])
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
        plt.savefig(self.output_file + 'isomer.png', format='png', dpi=300,bbox_inches='tight')
        plt.show()     
        
        self.write_file_nested_list(self.output_file + 'isomer.txt',isomer_list)
               
    def Read_NIST(self):
        """
        NIST.txtからデータを読み取る関数
        """
        self.nist_list = []
        self.natural_element_nist_list = []#自然界に存在する核種の元素記号
        self.natural_comp_list = []#自然界に存在する核種の存在比
        self.natural_ramass_list = []#自然界に存在する核種の相対原子質量
        self.natural_massnum_list = []#自然界に存在する核種の原子番号
        self.natural_saw_list = []#自然界に存在する核種の平均質量数
        self.element_nist_list = []#すべての元素記号
        self.massnum_list = []#すべての原子番号
        self.nist_list = self.Read_n_splited(self.input_file + self.NIST_file)

        for i in range(math.floor(len(self.nist_list)/8)):
            if len(self.nist_list[8*i+4]) != 25:#Isotopic Compositionのデータがあるとき（存在比が有意）
                if self.nist_list[8*i+2][15] == " ":#質量数が１桁の時
                    self.natural_element_nist_list.append(self.nist_list[8*i+1][16]+self.nist_list[8*i+1][17]+"-"+self.nist_list[8*i+2][14].rjust(3))
                elif self.nist_list[8*i+2][16] == " ":#質量数が2桁の時
                    self.natural_element_nist_list.append(self.nist_list[8*i+1][16]+self.nist_list[8*i+1][17]+"-"+(self.nist_list[8*i+2][14]+self.nist_list[8*i+2][15]).rjust(3))
                else:#質量数が３桁の時
                    self.natural_element_nist_list.append(self.nist_list[8*i+1][16]+self.nist_list[8*i+1][17]+"-"+self.nist_list[8*i+2][14]+self.nist_list[8*i+2][15]+self.nist_list[8*i+2][16])

                #同位体存在比を抽出--------------------
                s=self.nist_list[8*i+4].split()
                comp = ""
                for k in range(len(s[3])):
                    if s[3][k] == "(":
                        break
                    comp = comp + s[3][k]
                self.natural_comp_list.append(comp.ljust(10))

                #相対原子質量を抽出--------------------
                s = self.nist_list[8*i+3].split()
                ramass = ""
                for l in range(len(s[4])):
                    if s[4][l] == "(":
                        break
                    ramass = ramass + s[4][l]
                self.natural_ramass_list.append(ramass)

                #原子番号を抽出---------------------
                self.natural_massnum_list.append(int(self.nist_list[8*i][16]+self.nist_list[8*i][17]+self.nist_list[8*i][18]))

                #平均原子量を抽出
                s = self.nist_list[8*i+5].split()
                saw = ""
                for j in range(len(s[4])):
                    if s[4][j] == "(":
                        break
                    if s[4][j] != "[" and s[4][j] != "]":
                        saw = saw + s[4][j]
                if len(saw.split(",")) == 2:
                    d = saw.split(",")
                    saw = (float(d[0]) + float(d[1]))/2
                self.natural_saw_list.append(saw)


        #すべての元素を抽出------------------------------------------------------------------------------
        for i in range(math.floor(len(self.nist_list)/8)):
            self.element_nist_list.append(self.nist_list[8*i+1][16]+self.nist_list[8*i+1][17])

            #原子番号を抽出---------------------
            self.massnum_list.append(int(self.nist_list[8*i][16]+self.nist_list[8*i][17]+self.nist_list[8*i][18]))    
            
    def Plot_activation_cross_section(self):
        activation_cross_section_x = []
        activation_cross_section_y = []
        label_list = []
        self.Read_NIST()
        files = os.listdir(self.input_file+self.JENDL_file)
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
            num_lines = len(open(self.input_file+self.JENDL_file+t).readlines())#最終行を取得
            lines = self.Read_n_splited(self.input_file + self.JENDL_file + t)
            
            n=0
            for line in lines:
                n = n + 1
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
                MF10Q = self.numeric(text_list[k][11]+text_list[k][12]+text_list[k][13]+text_list[k][14]+text_list[k][15]+text_list[k][16]\
                                +text_list[k][17]+text_list[k][18]+text_list[k][19]+text_list[k][20]+text_list[k][21])
                
                #グランドステイトのQ値
                MF10Q0 = self.numeric(text_list[k][0]+text_list[k][1]+text_list[k][2]+text_list[k][3]+text_list[k][4]+text_list[k][5]\
                                 +text_list[k][6]+text_list[k][7]+text_list[k][8]+text_list[k][9]+text_list[k][10])
                
                #グランドステイトからの差分
                MF10dQ = (float(MF10Q0) - float(MF10Q))/1000
                data_list.append(MF10dQ)
                
                #Q値を追加-------------------------------------------------
                if float(MF10Q) < 0:#Q値が負の時（＝閾値がある）
                    #閾値
                    MF10Threshold = self.numeric(text_list[k+2][1]+text_list[k+2][2]+text_list[k+2][3]+text_list[k+2][4]+text_list[k+2][5]+text_list[k+2][6]\
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
                                        cross_data_list.append(self.numeric(text_list[u].split()[a][0]+text_list[u].split()[a][1]+text_list[u].split()[a][2]\
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
                                        cross_data_list.append(self.numeric(text_list[u].split()[a][0]+text_list[u].split()[a][1]+text_list[u].split()[a][2]+text_list[u].split()[a][3]\
                                                                       +text_list[u].split()[a][4]+text_list[u].split()[a][5]+text_list[u].split()[a][6]+text_list[u].split()[a][7]\
                                                                       +text_list[u].split()[a][8]+text_list[u].split()[a][9]))  
                cross_section_per_element_list.append(cross_data_list)#反応ごとの断面積スペクトルを核種でまとめる
                reaction_per_element_list.append(data_list)#それぞれの反応の情報を核種でまとめる
                
               
            
            if len(reaction_per_element_list) != 0:
                for i in range(len(reaction_per_element_list)):
                            
                    #自分が指定した反応の時
                    for iii in plot_activation_cross_section:
                        if reaction_per_element_list[i][0] == iii[0] and reaction_per_element_list[i][1] == iii[1] and reaction_per_element_list[i][2] == iii[2]:
                            #断面積データを抜き出す
                            x_list = []
                            y_list = []
                            for iiii in range(int(len(cross_section_per_element_list[i])/2)):
                                x_list.append(float(cross_section_per_element_list[i][2*iiii]))
                                y_list.append(float(cross_section_per_element_list[i][2*iiii+1]))
                            activation_cross_section_x.append(x_list)
                            activation_cross_section_y.append(y_list)
                            label_list.append('$\mathrm{^{%s}%s (%sQ=%s}$'% (reaction_per_element_list[i][0].split('-')[1],reaction_per_element_list[i][0].split('-')[0],reaction_per_element_list[i][1].split("(")[1][:5],reaction_per_element_list[i][2])) 
    
                        
        for i in range(len(activation_cross_section_x)):
            plt.plot(activation_cross_section_x[i],activation_cross_section_y[i],label = label_list[i])
            print(activation_cross_section_y[i])
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
        plt.savefig(self.output_file + 'Activation cross section.png', format='png', dpi=300,bbox_inches='tight')
        plt.show()
    
    def Make_reaction_list(self):
        """
        JENDLからデータを読み取り、NISTの情報と合わせてreaction_listをテキストファイルで出力する関数
        """
        cross_section_x_list = []#ｘとｙに分けた断面積スペクトル
        cross_section_y_list = []#ｘとｙに分けた断面積スペクトル
                                    
        #ファイル名を取得し、ファイルごとに情報を抜き出す
        files = os.listdir(self.input_file+self.JENDL_file)
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
            num_lines = len(open(self.input_file+self.JENDL_file+t).readlines())#最終行を取得
            lines = self.Read_n_splited(self.input_file + self.JENDL_file + t)
            
            n=0
            for line in lines:
                n = n + 1
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
                MF10Q = self.numeric(text_list[k][11]+text_list[k][12]+text_list[k][13]+text_list[k][14]+text_list[k][15]+text_list[k][16]\
                                +text_list[k][17]+text_list[k][18]+text_list[k][19]+text_list[k][20]+text_list[k][21])
                
                #グランドステイトのQ値
                MF10Q0 = self.numeric(text_list[k][0]+text_list[k][1]+text_list[k][2]+text_list[k][3]+text_list[k][4]+text_list[k][5]\
                                 +text_list[k][6]+text_list[k][7]+text_list[k][8]+text_list[k][9]+text_list[k][10])
                
                #グランドステイトからの差分
                MF10dQ = (float(MF10Q0) - float(MF10Q))/1000
                data_list.append(MF10dQ)
                
                #Q値を追加-------------------------------------------------
                if float(MF10Q) < 0:#Q値が負の時（＝閾値がある）
                    #閾値
                    MF10Threshold = self.numeric(text_list[k+2][1]+text_list[k+2][2]+text_list[k+2][3]+text_list[k+2][4]+text_list[k+2][5]+text_list[k+2][6]\
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
                                        cross_data_list.append(self.numeric(text_list[u].split()[a][0]+text_list[u].split()[a][1]+text_list[u].split()[a][2]\
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
                                        cross_data_list.append(self.numeric(text_list[u].split()[a][0]+text_list[u].split()[a][1]+text_list[u].split()[a][2]+text_list[u].split()[a][3]\
                                                                       +text_list[u].split()[a][4]+text_list[u].split()[a][5]+text_list[u].split()[a][6]+text_list[u].split()[a][7]\
                                                                       +text_list[u].split()[a][8]+text_list[u].split()[a][9]))  
                cross_section_per_element_list.append(cross_data_list)#反応ごとの断面積スペクトルを核種でまとめる
                reaction_per_element_list.append(data_list)#それぞれの反応の情報を核種でまとめる
            
            if len(reaction_per_element_list) != 0:
                for i in range(len(reaction_per_element_list)):
                    for ii in range(len(self.natural_element_nist_list)):
                        if self.natural_element_nist_list[ii] == reaction_per_element_list[i][0]:#自然界に存在するとき
                            
                            #閾値の条件を満たすとき
                            if float(reaction_per_element_list[i][3]) <= self.Threshold_MAX and float(reaction_per_element_list[i][3]) >= self.Threshold_MIN:
                                
                                ##じぶんで書いた以下の反応のif分に漏れがないか調べるためのprint----------------------------------------
                                ext = 0
                                for iii in ll:
                                    if iii in reaction_per_element_list[i][1]:
                                        ext =1
                                if ext != 1:
                                    print(reaction_per_element_list[i],"この反応が漏れています")#じぶんで書いた以下の反応のif分に漏れがないか調べるためのprint  
                                    
                                if  "(n,He3)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-2:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-2)#原子番号が2減る
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)#質量数が2減る
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                        
                                if  "(n,t)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-1:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-1)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                
                                if  "(n,2na)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-2:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-2)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-5)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                
                                if  "(n,2p)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-2:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-2)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,d)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-1:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-1)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,2n)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii])
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,np)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-1:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-1)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-1)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,p)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-1:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-1)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5]))
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                
                                if  "(n,pa)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-3:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-3)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-4)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,na)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-2:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-2)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-4)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,n')" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii])
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5]))
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,a)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-2:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-2)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-3)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,nd)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-1:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-1)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,2np)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-1:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-1)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
        
                                if  "(n,3n)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii])
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-2)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
                                            break
                                            
                                if  "(n,nt)" in reaction_per_element_list[i][1]:
                                    for j in range(len(self.massnum_list)):
                                        if self.massnum_list[j] == self.natural_massnum_list[ii]-1:
                                            reaction_per_element_list[i].append(self.natural_massnum_list[ii]-1)
                                            reaction_per_element_list[i].append(int(reaction_per_element_list[i][0][3]+reaction_per_element_list[i][0][4]+reaction_per_element_list[i][0][5])-3)
                                            reaction_per_element_list[i].append(self.element_nist_list[j].rstrip())
                                            reaction_per_element_list[i].append(self.natural_comp_list[ii].rstrip())
                                            reaction_per_element_list[i].append(self.natural_saw_list[ii])
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
                                
                                    
                                #断面積スペクトルを使って∫σΦを計算-------------------------------------------------
                                reaction_rate = 0
                                for iii in range(len(self.spectrum_x_list)):
                                    if self.Calculate_interpolation(x_list,y_list,self.spectrum_x_list[iii]*10**6,2) > 0:
                                        reaction_rate = reaction_rate + self.Calculate_interpolation(x_list,y_list,self.spectrum_x_list[iii]*10**6,2)*self.spectrum_y_list[iii]*10**(-24)
                                    else:
                                        reaction_rate = reaction_rate + 0
                                
                                #reaction_per_element_list[i].append(reaction_rate)     
                                self.reaction_list.append(reaction_per_element_list[i])
                    
                    if reaction_per_element_list[i][0] == "Nb- 93" and reaction_per_element_list[i][1] == "MF=10 MT= 16 (n,2n) reaction  " and \
                    reaction_per_element_list[i][2] == 135.5:
                        #reaction_per_element_listが条件を満たすときその断面積スペクトルを取得-------------------
                        x_list = []
                        y_list = []
                        
                        for iii in range(int(len(cross_section_per_element_list[i])/2)):
                            x_list.append(float(cross_section_per_element_list[i][2*iii]))
                            y_list.append(float(cross_section_per_element_list[i][2*iii+1]))
            
                        #断面積スペクトルを保存-----------------------------------
                        self.write_file_non_nested_list(self.output_file + '/Phase2/cross_section_Nb_x.txt', x_list)
                        self.write_file_non_nested_list(self.output_file + '/Phase2/cross_section_Nb_y.txt', y_list)                        
                    
        self.write_file_nested_list(self.output_file + '/Phase2/reaction_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt',self.reaction_list)
        self.write_file_nested_list(self.output_file + '/Phase2/cross_section_x' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt',cross_section_x_list)
        self.write_file_nested_list(self.output_file + '/Phase2/cross_section_y' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt',cross_section_y_list)
        
    def Make_result_list(self):
        reaction_list = self.Read_tab_splited(self.output_file + '/Phase2/reaction_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt' )
        cross_section_x_list = self.Read_tab_splited(self.output_file + '/Phase2/cross_section_x' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt' )
        cross_section_y_list = self.Read_tab_splited(self.output_file + '/Phase2/cross_section_y' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt' )
        density_list = self.Read_comma_splited(self.input_file + self.density_file)
        Gammas_list = self.Read_comma_splited(self.input_file + self.Gammas_file)
        cross_section_per_gammas_x_list = []
        cross_section_per_gammas_y_list = []
        
        for i in tqdm(range(len(reaction_list))):
            k = 0 #jendlから読み取った放射化後の順位のデータが自分で作ったgammmas.csvにあるかどうかを判別するための変数
            
            #reaction_listに密度を追加-------------------------------------------------------------------------------------------------
            pp = 0#密度データがdensity.csvにあるかどうかを判別するための変数
            for ii in range(len(density_list)):
                if density_list[ii][0].ljust(2) == reaction_list[i][0][0]+reaction_list[i][0][1]:
                    reaction_list[i].append(density_list[ii][1])
                    pp = 1
            if pp == 0:
                print("density.csvに",reaction_list[i][0][0]+reaction_list[i][0][1],"を追加してください")
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
        
            #                 if float(Gammas_list[ii][3]) < float(efficiency_list[0][0]):
            #                     print("検出効率が測定値以下の近似曲線で外挿されています",reaction_list[i],Gammas_list[ii][3])
        
        
                            exsistence_attenuation_file = 0#吸収係数のファイルがあるかどうかをみる変数
                            Attenuation_files = os.listdir(self.input_file + self.attenuation_file)
                            for q in Attenuation_files:
                                if q == (reaction_list[i][0][0]+reaction_list[i][0][1]).rstrip()+".txt":
                                    exsistence_attenuation_file = 1
                            if exsistence_attenuation_file == 0:
                                print("Attenuationに",(reaction_list[i][0][0]+reaction_list[i][0][1]).rstrip(),"がありません。追加してください。")
                                a.append("Attenuationのデータがありません")
                            else:        
                                #質量吸収係数を読み取る、a2に入れる-----------------------------------------------------------------------
                                a2 = self.Read_tab_splited(self.input_file + self.attenuation_file +'/%s.txt'%(reaction_list[i][0][0]+reaction_list[i][0][1]).rstrip())
                                
                                #質量吸収係数を算出するためにガンマ線のエネルギーが質量吸収係数のデータのどの座標にあるかを調べる--------------
                                attenuation = self.Calculate_interpolation([r[0] for r in a2], [r[1] for r in a2], float(Gammas_list[ii][3])*10**(-3), 2)

                            a.append(attenuation)

                            efficiency = self.func(float(Gammas_list[ii][3]),*[float(x) for x in self.Read_comma_splited(output_file+'Curve_fit.txt')[0]])
                            a.append(efficiency)
        
                        elif float(Gammas_list[ii][3]) == 0:
                            a.append(Gammas_list[ii][3])
                            a.append(Gammas_list[ii][4])
                            a.append(Gammas_list[ii][5])
                            a.append(0)
                            a.append(0)
        
                        self.result_list.append(a)
                        cross_section_per_gammas_x_list.append(cross_section_x_list[i])
                        cross_section_per_gammas_y_list.append(cross_section_y_list[i])
        
            if k == 0:#jendlから読み取ったデータがgammmas.csvにない
        #         print("jendlから読み取ったデータがgammmas.csvにありません。追加してください",reaction_list[i][4],reaction_list[i][5],reaction_list[i][6],reaction_list[i][2],reaction_list[i][1])
                print(reaction_list[i][4],reaction_list[i][5],reaction_list[i][6],reaction_list[i][2],reaction_list[i][1])

        self.write_file_nested_list(self.output_file + '/Phase3/result_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt',self.result_list)
        self.write_file_nested_list(self.output_file + '/Phase3/cross_section_x' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt',cross_section_per_gammas_x_list)
        self.write_file_nested_list(self.output_file + '/Phase3/cross_section_y' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt',cross_section_per_gammas_y_list)
      
    def Calculate_count(self):
        #スペクトルファイルを読み取り、ｘとｙに分ける--------------------------------------
        spectrum_list = self.Read_tab_splited(self.input_file + self.spectrum_file)
        for iii in range(len(spectrum_list)):
            self.spectrum_x_list.append(float(spectrum_list[iii][0]))   
            self.spectrum_y_list.append(float(spectrum_list[iii][1]))
        
        #[元素記号-原子番号、反応の種類、順位、閾値、娘核の原子番号、娘核の質量数、娘核の元素記号、存在比、平均原子量、密度、ガンマ線のエネルギー、放出比、半減期、質量吸収係数、検出効率]
        result_list = self.Read_tab_splited(self.output_file + '/Phase3/result_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_x_list = self.Read_tab_splited(self.output_file + '/Phase3/cross_section_x' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_y_list = self.Read_tab_splited(self.output_file + '/Phase3/cross_section_y' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_x_Nb_list = self.Read_n_splited(self.output_file + '/Phase2/cross_section_Nb_x.txt')
        cross_section_y_Nb_list = self.Read_n_splited(self.output_file + '/Phase2/cross_section_Nb_y.txt')
        N = 6.022*10**23
        

        for l in l_list:
            count_list = []
            info_list = []
            count_sum_list = []
            info_sum_list = []
            
            for i in tqdm(range(len(result_list))):
                count_per_list = []
                ticks = []
                
                if float(result_list[i][3]) <= self.Threshold_MAX_calc and float(result_list[i][3]) >= self.Threshold_MIN_calc\
                and not '*' in result_list[i][11] and float(result_list[i][10]) != 0 and result_list[i][12] != '0':
                    
                    info_list.append(result_list[i])
                    a = float(result_list[i][8])#平均原子量
                    b = float(result_list[i][7])#存在比
                    
                    #∫σΦdE
                    reaction_rate = 0
                    for ii in range(len(self.spectrum_x_list)):
                        if self.Calculate_interpolation(cross_section_x_list[i],cross_section_y_list[i],self.spectrum_x_list[ii]*10**6,2) > 0:
                            reaction_rate = reaction_rate + self.Calculate_interpolation(cross_section_x_list[i],cross_section_y_list[i],self.spectrum_x_list[ii]*10**6,2)\
                                *self.spectrum_y_list[ii]*10**(-24)
                        else:
                            reaction_rate = reaction_rate + 0
                    c = reaction_rate
                    
                    #半減期
                    d = ''
                    for ii in range(len(result_list[i][12])-1):
                        d = d + result_list[i][12][ii]
                    d = float(d)
                    if result_list[i][12][len(result_list[i][12])-1] == "s":
                        d = d
                    elif result_list[i][12][len(result_list[i][12])-1] == "m":
                        d = d*60
                    elif result_list[i][12][len(result_list[i][12])-1] == "h":
                        d = d*60*60
                    elif result_list[i][12][len(result_list[i][12])-1] == "d":
                        d = d*60*60*24
                    elif result_list[i][12][len(result_list[i][12])-1] == "y":
                        d = d*60*60*24*365
                    d = float(d)
                    
                    e = ''
                    if '<' in result_list[i][11] or '~' in result_list[i][11] or '>' in result_list[i][11]:
                        for ii in range(1,len(result_list[i][11])):
                            e = e + result_list[i][11][ii]
                    else:
                        e = float(result_list[i][11])
                    e = float(e)/100
                    
                    f = float(result_list[i][14])#検出効率
                    g = float(result_list[i][13])#質量吸収係数
                    h = float(result_list[i][9])#平均密度
                    #print(e)
                    for ii in range(self.day):
                        if ii == 0:
                            housyaka =  b * (h*(r*r*3.14)*l)/a * N * c * d/np.log(2) * (1-np.exp(-1*np.log(2)/d*ta*3600))*DT
                            not_housyaka =  b * (h*(r*r*3.14)*l)/a * N - housyaka
                        else:
                            not_housyaka = not_housyaka - not_housyaka * c * d/np.log(2) * (1-np.exp(-1*np.log(2)/d*ta*3600))*DT
                            housyaka = not_housyaka * c * d/np.log(2) * (1-np.exp(-1*np.log(2)/d*ta*3600))*DT + housyaka * np.exp(-1*np.log(2)/d*(tb+ta)*3600)
                        for t in range(1,self.tc):
                            count =  housyaka * (1-np.exp(-1*l*h*g))/(l*h*g) * e * f  * np.exp(-1*np.log(2)/d*1*3600) * (1-np.exp(-1*np.log(2)/d*t*24*3600))
                            count_per_list.append(count)
                    count_list.append(count_per_list)

            
            self.write_file_nested_list(self.output_file + '/Calculate_count/count_length=' + str(l) + '.txt',count_list)
            self.write_file_nested_list(self.output_file + '/Calculate_count/info_length=' + str(l) + '.txt',info_list)
            
            info_duplicate = []
            for i in tqdm(range(len(info_list))):
                #info_duplicate.append(info_list[i][0] + ',' + info_list[i][1])
                a = []
                a.append(info_list[i][0])
                a.append(info_list[i][1])
                a.append(info_list[i][3])
                info_duplicate.append(a)
            #info_not_duplicate = list(set(info_duplicate))
            info_not_duplicate = self.get_unique_list(info_duplicate)

            for i in tqdm(range(len(info_not_duplicate))):
                count_sum = [0]*len(count_list[1])
                info_sum = info_not_duplicate[i]
                for ii in range(len(info_list)):
                    if info_list[ii][0] == info_not_duplicate[i][0] and info_list[ii][1] == info_not_duplicate[i][1] and info_list[ii][3] == info_not_duplicate[i][2]:
                        count_sum = [x + y for (x,y) in zip(count_sum,count_list[ii])]
                        info_sum.append(info_list[ii][2]+','+info_list[ii][10]+','+info_list[ii][11])
                count_sum_list.append(count_sum)
                info_sum_list.append(info_sum)
            self.write_file_nested_list(self.output_file+'/Calculate_count/count_sum_length=' + str(l) + '.txt',count_sum_list)
            self.write_file_nested_list(self.output_file+'/Calculate_count/info_sum_length=' + str(l) + '.txt',info_sum_list)
                
            #Nbのカウント数の計算
            Nb_a = 92.90637
            Nb_b = 1.0
            Nb_d = 10.15*24*60*60
            Nb_e = 0.9907
            Nb_f = float("5.61E-03")
            Nb_g = 0.0609
            Nb_h = 8.57
            count_Nb_list = []
            reaction_rate = 0
            
            for ii in range(len(self.spectrum_x_list)):
                if self.Calculate_interpolation(cross_section_x_Nb_list,cross_section_y_Nb_list,self.spectrum_x_list[ii]*10**6,2) > 0:
                    reaction_rate = reaction_rate + self.Calculate_interpolation(cross_section_x_Nb_list,cross_section_y_Nb_list,self.spectrum_x_list[ii]*10**6,2)*self.spectrum_y_list[ii]*10**(-24)
                else:
                    reaction_rate = reaction_rate + 0
            sigma_reaction_Nb = reaction_rate  
            
            for i in range(self.day):
                if i == 0:
                    housyaka =  Nb_b * (Nb_h*(r*r*3.14)*l)/Nb_a * N * sigma_reaction_Nb * Nb_d/np.log(2) * (1-np.exp(-1*np.log(2)/Nb_d*ta*3600))*DT
                    not_housyaka =  Nb_b * (Nb_h*(r*r*3.14)*l)/Nb_a * N - housyaka
                else:
                    not_housyaka = not_housyaka - not_housyaka * sigma_reaction_Nb * Nb_d/np.log(2) * (1-np.exp(-1*np.log(2)/Nb_d*ta*3600))*DT
                    housyaka = not_housyaka * sigma_reaction_Nb * Nb_d/np.log(2) * (1-np.exp(-1*np.log(2)/Nb_d*ta*3600))*DT + housyaka * np.exp(-1*np.log(2)/Nb_d*(tb+ta)*3600)
                for t in range(1,self.tc):
                    count_Nb =  housyaka * (1-np.exp(-1*l*Nb_h*Nb_g))/(l*Nb_h*Nb_g) * Nb_e * Nb_f  * np.exp(-1*np.log(2)/Nb_d*1*3600) * (1-np.exp(-1*np.log(2)/Nb_d*t*24*3600))
                    count_Nb_list.append(count_Nb)
            
            self.write_file_non_nested_list(self.output_file + '/Calculate_count/count_Nb_length=' + str(l) + '.txt', count_Nb_list)         
      
    def Calculate_count_self_absorbtion(self):
        #スペクトルファイルを読み取り、ｘとｙに分ける--------------------------------------
        spectrum_list = self.Read_tab_splited(self.input_file + self.spectrum_file)
        for iii in range(len(spectrum_list)):
            self.spectrum_x_list.append(float(spectrum_list[iii][0]))   
            self.spectrum_y_list.append(float(spectrum_list[iii][1]))
        
        #[元素記号-原子番号、反応の種類、順位、閾値、娘核の原子番号、娘核の質量数、娘核の元素記号、存在比、平均原子量、密度、ガンマ線のエネルギー、放出比、半減期、質量吸収係数、検出効率]
        result_list = self.Read_tab_splited(self.output_file + '/Phase3/result_list_' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_x_list = self.Read_tab_splited(self.output_file + '/Phase3/cross_section_x' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_y_list = self.Read_tab_splited(self.output_file + '/Phase3/cross_section_y' + str(self.Threshold_MIN) + '_' + str(self.Threshold_MAX) + '.txt')
        cross_section_x_Nb_list = self.Read_n_splited(self.output_file + '/Phase2/cross_section_Nb_x.txt')
        cross_section_y_Nb_list = self.Read_n_splited(self.output_file + '/Phase2/cross_section_Nb_y.txt')
        N = 6.022*10**23
        

        for l in l_list:
            count_list = []
            info_list = []
            count_sum_list = []
            info_sum_list = []
            
            for i in tqdm(range(len(result_list))):
                count_per_list = []
                ticks = []
                

                if result_list[i][0] in calc_list\
                    and not '*' in result_list[i][11] and float(result_list[i][10]) != 0 and result_list[i][12] != '0':
                
                    info_list.append(result_list[i])
                    
                    #全断面積
                    ff = open(self.input_file + 'Considering_neutron_absorbtion/' + result_list[i][0].split("-")[0]+result_list[i][0].split("-")[1]+".txt")
                    line = ff.readline()
                    all_cross_section_x_list = []
                    all_cross_section_y_list = []
                    p = 0
                    while line:
                        p = p + 1
                        if p > 4:
                            for ii in range(int((len(line[:-15].split()))/2)):
                                all_cross_section_x_list.append(float(self.numeric(line[:-15].split()[ii*2])))
                                all_cross_section_y_list.append(float(self.numeric(line[:-15].split()[ii*2+1])))
                        line = ff.readline()
                    ff.close
                    
                    a = float(result_list[i][8])#平均原子量
                    b = float(result_list[i][7])#存在比
                    
                    
                    #半減期
                    d = ''
                    for ii in range(len(result_list[i][12])-1):
                        d = d + result_list[i][12][ii]
                    d = float(d)
                    if result_list[i][12][len(result_list[i][12])-1] == "s":
                        d = d
                    elif result_list[i][12][len(result_list[i][12])-1] == "m":
                        d = d*60
                    elif result_list[i][12][len(result_list[i][12])-1] == "h":
                        d = d*60*60
                    elif result_list[i][12][len(result_list[i][12])-1] == "d":
                        d = d*60*60*24
                    elif result_list[i][12][len(result_list[i][12])-1] == "y":
                        d = d*60*60*24*365
                    d = float(d)
                    
                    e = ''
                    if '<' in result_list[i][11] or '~' in result_list[i][11] or '>' in result_list[i][11]:
                        for ii in range(1,len(result_list[i][11])):
                            e = e + result_list[i][11][ii]
                    else:
                        e = float(result_list[i][11])
                    e = float(e)/100
                    
                    f = float(result_list[i][14])#検出効率
                    g = float(result_list[i][13])#質量吸収係数
                    h = float(result_list[i][9])#平均密度
                    
                    #∫σΦdE
                    reaction_rate = 0
                    for ii in range(len(self.spectrum_x_list)):
                        if self.Calculate_interpolation(cross_section_x_list[i],cross_section_y_list[i],self.spectrum_x_list[ii]*10**6,2) > 0:
                            reaction_rate = reaction_rate + self.Calculate_interpolation(cross_section_x_list[i],cross_section_y_list[i],self.spectrum_x_list[ii]*10**6,2)\
                                *self.spectrum_y_list[ii]*10**(-24)*(1-np.exp(-h/a*6*10**23*\
                                self.Calculate_interpolation(all_cross_section_x_list,all_cross_section_y_list,self.spectrum_x_list[ii]*10**6,2)*10**(-24)*l))/\
                                (h/a*6*10**23*self.Calculate_interpolation(all_cross_section_x_list,all_cross_section_y_list,self.spectrum_x_list[ii]*10**6,2)*10**(-24)*l)

                        else:
                            reaction_rate = reaction_rate + 0
                    c = reaction_rate

                    
                    
                    for ii in range(self.day):
                        if ii == 0:
                            housyaka =  b * (h*(r*r*3.14)*l)/a * N * c * d/np.log(2) * (1-np.exp(-1*np.log(2)/d*ta*3600))*DT
                            not_housyaka =  b * (h*(r*r*3.14)*l)/a * N - housyaka
                        else:
                            not_housyaka = not_housyaka - not_housyaka * c * d/np.log(2) * (1-np.exp(-1*np.log(2)/d*ta*3600))*DT
                            housyaka = not_housyaka * c * d/np.log(2) * (1-np.exp(-1*np.log(2)/d*ta*3600))*DT + housyaka * np.exp(-1*np.log(2)/d*(tb+ta)*3600)
                        for t in range(1,self.tc):
                            count =  housyaka * (1-np.exp(-1*l*h*g))/(l*h*g) * e * f  * np.exp(-1*np.log(2)/d*1*3600) * (1-np.exp(-1*np.log(2)/d*t*24*3600))
                            count_per_list.append(count)
                    count_list.append(count_per_list)

            
            self.write_file_nested_list(self.output_file + '/Calculate_count_self_absorbtion/count_length=' + str(l) + '.txt',count_list)
            self.write_file_nested_list(self.output_file + '/Calculate_count_self_absorbtion/info_length=' + str(l) + '.txt',info_list)
            
            info_duplicate = []
            for i in tqdm(range(len(info_list))):
                #info_duplicate.append(info_list[i][0] + ',' + info_list[i][1])
                a = []
                a.append(info_list[i][0])
                a.append(info_list[i][1])
                a.append(info_list[i][3])
                info_duplicate.append(a)
            #info_not_duplicate = list(set(info_duplicate))
            info_not_duplicate = self.get_unique_list(info_duplicate)

            for i in tqdm(range(len(info_not_duplicate))):
                count_sum = [0]*len(count_list[1])
                info_sum = info_not_duplicate[i]
                for ii in range(len(info_list)):
                    if info_list[ii][0] == info_not_duplicate[i][0] and info_list[ii][1] == info_not_duplicate[i][1] and info_list[ii][3] == info_not_duplicate[i][2]:
                        count_sum = [x + y for (x,y) in zip(count_sum,count_list[ii])]
                        info_sum.append(info_list[ii][2]+','+info_list[ii][10]+','+info_list[ii][11])
                count_sum_list.append(count_sum)
                info_sum_list.append(info_sum)
            self.write_file_nested_list(self.output_file+'/Calculate_count_self_absorbtion/count_sum_length=' + str(l) + '.txt',count_sum_list)
            self.write_file_nested_list(self.output_file+'/Calculate_count_self_absorbtion/info_sum_length=' + str(l) + '.txt',info_sum_list)
                
            #Nbのカウント数の計算
            Nb_a = 92.90637
            Nb_b = 1.0
            Nb_d = 10.15*24*60*60
            Nb_e = 0.9907
            Nb_f = float("5.61E-03")
            Nb_g = 0.0609
            Nb_h = 8.57
            count_Nb_list = []
            reaction_rate = 0
            
            for ii in range(len(self.spectrum_x_list)):
                if self.Calculate_interpolation(cross_section_x_Nb_list,cross_section_y_Nb_list,self.spectrum_x_list[ii]*10**6,2) > 0:
                    reaction_rate = reaction_rate + self.Calculate_interpolation(cross_section_x_Nb_list,cross_section_y_Nb_list,self.spectrum_x_list[ii]*10**6,2)*self.spectrum_y_list[ii]*10**(-24)
                else:
                    reaction_rate = reaction_rate + 0
            sigma_reaction_Nb = reaction_rate  
            
            for i in range(self.day):
                if i == 0:
                    housyaka =  Nb_b * (Nb_h*(r*r*3.14)*l)/Nb_a * N * sigma_reaction_Nb * Nb_d/np.log(2) * (1-np.exp(-1*np.log(2)/Nb_d*ta*3600))*DT
                    not_housyaka =  Nb_b * (Nb_h*(r*r*3.14)*l)/Nb_a * N - housyaka
                else:
                    not_housyaka = not_housyaka - not_housyaka * sigma_reaction_Nb * Nb_d/np.log(2) * (1-np.exp(-1*np.log(2)/Nb_d*ta*3600))*DT
                    housyaka = not_housyaka * sigma_reaction_Nb * Nb_d/np.log(2) * (1-np.exp(-1*np.log(2)/Nb_d*ta*3600))*DT + housyaka * np.exp(-1*np.log(2)/Nb_d*(tb+ta)*3600)
                for t in range(1,self.tc):
                    count_Nb =  housyaka * (1-np.exp(-1*l*Nb_h*Nb_g))/(l*Nb_h*Nb_g) * Nb_e * Nb_f  * np.exp(-1*np.log(2)/Nb_d*1*3600) * (1-np.exp(-1*np.log(2)/Nb_d*t*24*3600))
                    count_Nb_list.append(count_Nb)
            
            self.write_file_non_nested_list(self.output_file + '/Calculate_count_self_absorbtion/count_Nb_length=' + str(l) + '.txt', count_Nb_list)     
          
    def plot_count(self,folder):
        markers1 = [",","o","v","^","<",">"]
        for l in l_list:
            info_list = self.Read_tab_splited(self.output_file + str(folder) +'/info_sum_length=' + str(l) + '.txt')
            count_list = self.Read_tab_splited(self.output_file + str(folder) + '/count_sum_length=' + str(l) + '.txt')
            count_Nb_list = self.Read_n_splited(self.output_file + str(folder) + '/count_Nb_length=' + str(l) + '.txt')
            plot_list = []
            plot_info_list = []

            
            for i in range(len(count_list)):
                for ii in range(len(count_list[i])):
                    count_list[i][ii] = float(count_list[i][ii])
                   
            q = 0
            ticks = []
            for ii in range(self.day):
                for iii in range(1,self.tc):
                    ticks.append(q)
                    q = q + 1
            
            for i in tqdm(range(len(info_list)-1)):
                if float(info_list[i][2]) <= self.Threshold_MAX_plot and float(info_list[i][2]) >= self.Threshold_MIN_plot:
                    over = 0
                    for ii in range(len(count_list[i])):
                        if count_list[i][ii] > self.plot_lower_limit:
                            over = 1
                    if over == 1:
                        plt.scatter(ticks,count_list[i],label=info_list[i],s = 15,marker=markers1[i % 6])
                        plot_list.append(count_list[i])
                        plot_info_list.append(info_list[i])
    
                   
            plt.yscale('log')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
            plt.ylim(100,30000)
            plt.savefig(self.output_file + '/plot/The number of counts l =' + str(l) + '.png', format='png', dpi=300,bbox_inches='tight')
            plt.show()
            
            self.write_file_nested_list(self.output_file + '/plot/plot.txt',plot_list)
            self.write_file_nested_list(self.output_file + '/plot/info.txt',plot_info_list)
        

                
                    
if __name__ == "__main__":
    #Pathを指定
    input_file = './reference/'
    output_file = './OUTPUT/'
    
    
    #Phase_1
    efficiency_file = 'efficiency.txt'
    
    #Phase_2
    Threshold_MIN = 0E+6
    Threshold_MAX = 14.1E+6
    NIST_file = 'NIST.txt'
    JENDL_file = 'jendl-ad2017_300K.tar/jendl-ad2017_300K/'
    
    #Phase_3
    density_file = 'density.csv'
    Gammas_file = 'Gammas_molded.csv'
    attenuation_file = 'Attenuation'
    
    #Calculate_count
    spectrum_file = 'Fe_spectrum.txt'
    r = 1.5 #[cm] 箔の半径
    DT = 5*10**9
    ta = 8
    tb = 1
    l_list = [0.5] #[cm]箔の厚み
    day = 1 #何日間照射するか
    tc = 16 #計測時間
    Threshold_MIN_calc = 0
    Threshold_MAX_calc = 20E+6
    
    #Calculate_count_self_absorbtion
    calc_list = ['Au-197','Ir-191','Ir-193','Os-192','Re-187','Sb-123','Ta-181','Nb- 93']
    
    #Plot_count
    Threshold_MIN_plot = 0
    Threshold_MAX_plot = 20E+6
    plot_lower_limit = 1000
    
    
    #Plot_activation_cross_section
    #plot_activation_cross_section = [['Al- 27', 'MF=10 MT=107 (n,a) reaction   ',0.0],['Nb- 93','MF=10 MT= 16 (n,2n) reaction  ',135.5]]
    #plot_activation_cross_section = [['Al- 27', 'MF=10 MT=107 (n,a) reaction   ',0.0]]
    #plot_activation_cross_section = [['Yb-170', 'MF=10 MT= 16 (n,2n) reaction  ',0.0],['Yb-170', 'MF=10 MT= 16 (n,2n) reaction  ',24.2]]
    #plot_activation_cross_section = [['Hg-196', 'MF=10 MT= 16 (n,2n) reaction  ',0.0],['Hg-196', 'MF=10 MT= 16 (n,2n) reaction  ',176.08]]
    #plot_activation_cross_section = [['Os-192', 'MF=10 MT= 16 (n,2n) reaction  ',0.0],['Os-192', 'MF=10 MT= 16 (n,2n) reaction  ',74.38]]
    #plot_activation_cross_section = [['Xe-128', 'MF=10 MT= 16 (n,2n) reaction  ',0.0],['Xe-128', 'MF=10 MT= 16 (n,2n) reaction  ',297.1]]
    plot_activation_cross_section = [['In-115', "MF=3 MT= 7 (n,n') reaction  ",0.0]]
    #plot_activation_cross_section = [['Hg-196', 'MF=10 MT= 16 (n,2n) reaction  ',0.0],['Hg-196', 'MF=10 MT= 16 (n,2n) reaction  ',176.08]]
    #plot_activation_cross_section = [['Kr- 80', 'MF=10 MT= 16 Partial (n,2n) re',0.0],['Kr- 80', 'MF=10 MT= 16 Partial (n,2n) re',129.78]]
    
    
        
    execute = Optimum_foil(input_file,output_file,efficiency_file,NIST_file,JENDL_file,density_file,Gammas_file,attenuation_file,spectrum_file,\
                           Threshold_MIN,Threshold_MAX,Threshold_MIN_calc,Threshold_MAX_calc,Threshold_MIN_plot,Threshold_MAX_plot,day,tc,\
                               plot_activation_cross_section,plot_lower_limit)

    execute.Plot_activation_cross_section()
    #execute.Calculate_count()
    #execute.Calculate_count_self_absorbtion()
    #execute.plot_count('Calculate_count')
    #execute.plot_isomer() 
    #execute.Phase_3()



    
    

