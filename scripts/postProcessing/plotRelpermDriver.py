from matplotlib import pyplot as plt
from scipy import interpolate
import numpy as np
import argparse

def populate_input(fname):
    T = np.loadtxt(fname)
    return T

def ternatry_figure(int):
    h = plt.figure()
    h.clf()
    plt.axis('off')
    plt.plot([0,1/np.sqrt(3), 2/np.sqrt(3),0],[0,1,0,0],'-k')
    plt.axis('square')

    # add minor lines ?
    plt.text(1/np.sqrt(3)-0.03, 1+0.03, '$S_w$')
    plt.text(2/np.sqrt(3)+0.01, 0-0.03, '$S_o$')
    plt.text(0-0.05, 0-0.03, '$S_g$')

    h.tight_layout()
    return h

def main(fname):

    T = populate_input(fname[0])
    isThreePhase = False
    if T.shape[1] == 5 : # two phase
        Sw = T[:,3]#X
        Sg = T[:,2]#Y
    elif T.shape[1] == 7 :
        isThreePhase = True
        Sw = T[:,2]#X
        Sg = T[:,1]#Y
    else:
        raise NotImplemented()

    if(isThreePhase):

        #some preporcessing
        # Kro = interpolate.interp2d( Sw, Sg, T[:,4], kind='linear' )
        # X,Y = np.meshgrid(Sw,Sg)
        # X_ternary = ( (1-Y[:]) / np.cos(np.pi/6) ) - ( X[:] / np.tan(np.pi/3) )
        # X_ternary = np.reshape(X_ternary,(X.shape[0],X.shape[1]))
        # Y_ternary = X

    #try #2
        npt = int(np.sqrt(T[:,4].shape[0]))
        print(f'npt {npt}, shape{T[:,4].shape[0]}')
        Kro_ = np.reshape(T[:,4],(npt,npt))
        Sw = np.reshape(Sw,(npt,npt))
        Sg = np.reshape(Sg,(npt,npt))
        X_ternary = ( (1-Sg[:]) / np.cos(np.pi/6) ) - ( Sw[:] / np.tan(np.pi/3) )
        X_ternary = np.reshape(X_ternary,(npt,npt))
        Y_ternary = Sw


        h = ternatry_figure(1)
        kro_ = np.empty((X_ternary.shape[0],X_ternary.shape[1],))
        kro_[:] = np.nan
        for i in range(X_ternary.shape[0]):
            for j in range(X_ternary.shape[1]):
                if(1-Sw[i,j]-Sg[i,j]>=0):
                    # kro_[i,j] = Kro(X_ternary[i,j], Y_ternary[i,j])
                    kro_[i,j] = Kro_[i,j]

        plt.contourf(X_ternary,Y_ternary,kro_,levels=np.asarray([0,0.001,0.01,0.05,0.1,0.2,0.3,0.5,0.6,0.7,0.8,0.9,1.0]))
        # plt.show()
        plt.title('kro')
        plt.tight_layout()
        plt.savefig("relpermDriver_kro.pdf")

        plt.subplot(1,2,1)
        plt.plot(T[:-2,3],T[:-2,6],'or')
        plt.legend(['krw'])
        plt.xlabel('Sw')
        plt.ylabel('krw')
        plt.subplot(1,2,2)
        plt.plot(T[:-2,2],T[:-2,5],'ok')
        plt.legend(['krnw'])
        plt.xlabel('Snw')
        plt.ylabel('krnw')
        # plt.show()
        plt.tight_layout()
        plt.savefig("relpermDriver_kwg.pdf")
    else:
        plt.subplot(1,2,1)
        plt.plot(T[:-2,1],T[:-2,3],'or')
        plt.legend(['krw'])
        plt.xlabel('Sw')
        plt.ylabel('krw')
        plt.subplot(1,2,2)
        plt.plot(T[:-2,2],T[:-2,4],'ok')
        plt.legend(['krnw'])
        plt.xlabel('Snw')
        plt.ylabel('krnw')
        # plt.show()
        plt.tight_layout()
        plt.savefig("relpermDriver_kwg.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', type=str, nargs='+')
    args = parser.parse_args()
    main(args.files)




