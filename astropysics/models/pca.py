#Copyright (c) 2008 Erik Tollerud (etolleru@uci.edu)

class Pca(object):
    """
    A basic class for Principal Component Analysis (PCA).
    
    p is the number of dimensions, while N is the number of data points
    """
    _colors=('r','g','b','c','y','m','k') #defaults
    
    def __calc(self):
        A = self.A
        M=A-np.mean(A,axis=0)
        N=M/np.std(M,axis=0)
        
        self.M = M
        self.N = N
        self._eig = None
    
    def __init__(self,data,names=None):
        """
        p X N matrix input
        """
        from warnings import warn
        A = np.array(data).T
        n,p = A.shape
        self.n,self.p = n,p
        if p > n:
            warn('p > n - intentional?')
        self.A = A
        self._origA=A.copy()
            
        self.__calc()
        
        self._colors= np.tile(self._colors,int((p-1)/len(self._colors))+1)[:p]
        if names is not None and len(names) != p:
            raise ValueError('names must match data dimension')
        self.names = None if names is None else tuple([str(n) for n in names])
        
        
    def cov(self):
        return np.cov(self.N.T)
        
    def eig(self):
        if self._eig is None:
            res = np.linalg.eig(self.cov())
            sorti=np.argsort(res[0])[::-1]
            res=(res[0][sorti],res[1][:,sorti])
            self._eig=res
        return self._eig
    def evals(self):
        return self.eig()[0]
    def evecs(self):
        return self.eig()[1]
    
    def energies(self):
        v=self.evals()
        return v/np.sum(v)
        
    def plot2d(self,ix=0,iy=1,clf=True):
        import matplotlib.pyplot as plt
        x,y=self.N[:,ix],self.N[:,iy]
        if clf:
            plt.clf()
        plt.scatter(x,y)
        vals,evs=self.eig()
        #evx,evy=evs[:,ix],evs[:,iy]
        xl,xu=plt.xlim()
        yl,yu=plt.ylim()
        dx,dy=(xu-xl),(yu-yl)
        for val,vec,c in zip(vals,evs.T,self._colors):
            plt.arrow(0,0,val*vec[ix],val*vec[iy],head_width=0.05*(dx*dy/4)**0.5,fc=c,ec=c)
        #plt.arrow(0,0,vals[ix]*evs[ix,ix],vals[ix]*evs[iy,ix],head_width=0.05*(dx*dy/4)**0.5,fc='g',ec='g')
        #plt.arrow(0,0,vals[iy]*evs[ix,iy],vals[iy]*evs[iy,iy],head_width=0.05*(dx*dy/4)**0.5,fc='r',ec='r')
        if self.names is not None:
            plt.xlabel('$'+self.names[ix]+'/\\sigma$')
            plt.ylabel('$'+self.names[iy]+'/\\sigma$')
        
    def plot3d(self,ix=0,iy=1,iz=2,clf=True):
        import enthought.mayavi.mlab as M
        if clf:
            M.clf()
        z3=np.zeros(3)
        v=(self.evecs()*self.evals())
        M.quiver3d(z3,z3,z3,v[ix],v[iy],v[iz],scale_factor=5)
        M.points3d(self.N[:,ix],self.N[:,iy],self.N[:,iz],scale_factor=0.3)
        if self.names:
            M.axes(xlabel=self.names[ix]+'/sigma',ylabel=self.names[iy]+'/sigma',zlabel=self.names[iz]+'/sigma')
        else:
            M.axes()
        
    def sigclip(self,sigs):
        if np.isscalar(sigs):
            sigs=sigs*np.ones(self.N.shape[1])
        n = self.N.shape[0]
        m = np.all(np.abs(self.N) < sigs,axis=1)
        self.A=self.A[m]
        self.__calc()
        return n-sum(m)
        
    def reset(self):
        self.A = self._origA.copy()
        self.__calc()
        
        
    def project(self,enthresh=None,nPCs=None,cumen=None,vals=None):
        """
        projects the normalized values onto the components
        
        enthresh, nPCs, and cumen determine how many PCs to use
        
        if vals is None, the normalized data vectors are the values to project
        
        returns n,p(>threshold) dimension array
        """
        nonnones = sum([e != None for e in (enthresh,nPCs,cumen)])
        if nonnones == 0:
            m = slice(None)
        elif nonnones > 1:
            raise ValueError("can't specify more than one threshold")
        else:
            if enthresh is not None:
                m = self.energies() > enthresh
            elif nPCs is not None:
                m = slice(None,nPCs)
            elif cumen is not None:
                m = np.cumsum(self.energies()) <  cumen
            else:
                raise RuntimeError('Should be unreachable')
        
        if vals is None:
            vals = self.N.T
        else:
            if self.N.T.shape[0] != vals.shape[0]:
                raise ValueError("shape for vals doesn't match")
        proj = np.matrix(self.evecs()).T*vals
        return proj[m].T
            
    def deproject(self,A,normed=True):
        """
        input is an n X q array, where q <= p
        
        output is p X n
        """
        A=np.atleast_2d(A)
        n,q = A.shape
        p = self.A.shape[1]
        if q > p :
            raise ValueError("q > p")
        
        evinv=np.linalg.inv(np.matrix(self.evecs()).T)
        
        zs = np.zeros((n,p))
        zs[:,:q]=A
        
        proj = evinv*zs.T
        
        if normed:
            return np.array(proj.T).T
        else:
            mns=np.mean(self.A,axis=0)
            sds=np.std(self.M,axis=0)
            return (np.array(proj.T)*sds+mns).T
    
    def subtractPC(self,pc,vals=None):
        """
        pc can be a scalar or any sequence of pc indecies
        
        if vals is None, the source data is self.A, else whatever is in vals
        (which must be p x m)
        """
        if vals is None:
            vals = self.A
        else:
            vals = vals.T
            if vals.shape[1]!= self.A.shape[1]:
                raise ValueError("vals don't have the correct number of components")
        
        pcs=self.project()
        zpcs=np.zeros_like(pcs)
        zpcs[:,pc]=pcs[:,pc]
        upc=self.deproject(zpcs,False)
        
        A = vals.T-upc
        B = A.T*np.std(self.M,axis=0)
        return B+np.mean(self.A,axis=0)

