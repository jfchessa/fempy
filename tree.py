class tree(object):
    
    def __init__(self,id):
        self.id=id
        self.children = []
        self.parents = []
        
    def __repr__(self):
        return self.lrepr()
 
    def lrepr(self,l=0):
        ofs=l*'   '+'|--> '
        r = ofs+'node '+str(self.id)+' with '+str(len(self.children))+' children\n'
        for n in self.children:
            r = r + n.lrepr(l+1)
        return r    
            
                      
    def numparents(self):
        return len(self.parents)
        
    def numchildren(self):
        return len(self.children)
        
    def root(self):
        if self.numparents() == 0:
            return True
        else:
            return False
        
    def leaf(self):
        if self.numchildren() == 0:
            return True
        else:
            return False     

    def addbranch(self,node):
        self.children.append(node)
        node.parents.append(self)
        
    def numdecendents(self):
        c = len(self.children)
        for n in self.children:
            c = c + n.numdecendents()
        return c       
        
n0 = tree(0)  
n1 = tree(1)  
n2 = tree(2)

n0.addbranch(n1)
n0.addbranch(n2)

n3 = tree(3)
n2.addbranch(n3)

n3.addbranch(n1)

print n0

