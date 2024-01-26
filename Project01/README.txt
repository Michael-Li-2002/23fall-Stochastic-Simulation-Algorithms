ASA Project 1 李佳 2100010793
code:
  case_2d: 二维情形
    主程序: a_EnergyHeat_c_Correlation:  计算(a)(c)中的内能、比热、相关函数
                b_Magnetization:                     计算(b)中的磁化强度
                c_d_CorrelationLength:           计算(c)(d)中的关联长度
    函数:    metropolis:                               传统Metropolis方法, Single Spin Proposal
                Wolff:                                       Wolff算法(聚类)
                Hamilton:                                 计算Hamiltonian
                Correlation:                              计算\sigma_i \sigma_j, j-i=(\pm k,0), (0,\pm k)
                LeastSquare:                            最小二乘作线性拟合
  case_3d: 三维情形
    主程序: a_EnergyHeat_c_Correlation_3d:  计算(a)(c)中的内能、比热、相关函数
                b_Magnetization_3d:               计算(b)中的磁化强度
                c_d_CorrelationLength_3d:      计算(c)(d)中的关联长度
    函数:    metropolis:                               传统Metropolis方法, Single Spin Proposal
                Wolff:                                       Wolff算法(聚类)
                Hamilton:                                 计算Hamiltonian
                Correlation:                              计算\sigma_i \sigma_j, j-i=(\pm k,0,0), (0,\pm k,0),(0,0,\pm k)
                LeastSquare:                            最小二乘作线性拟合    