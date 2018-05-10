###################################################################
#
# Calculate exact pdf for an axis and allies land battle
#
# James Lischeske
#
# 5/8/18
#
####################################################################

import numpy as np
from scipy.misc import comb
import time


def gen_dice_pdf(N,diceLim):
    P = diceLim/6
    return np.asarray([comb(N,i)* P**i * (1-P)**(N-i) for i in range(N+1)])

def generalized_battle_pdf(nOnes, nTwos, nThrees, nFours):
    pdfOnes = gen_dice_pdf(nOnes, 1)
    pdfTwos = gen_dice_pdf(nTwos, 2)
    pdfThrees = gen_dice_pdf(nThrees, 3)
    pdfFours = gen_dice_pdf(nFours, 4)

    pdfFinal = np.zeros(len(pdfOnes)+len(pdfTwos)+len(pdfThrees)+len(pdfFours)-4+1)
    for i in range(len(pdfOnes)):
        for j in range(len(pdfTwos)):
            for k in range(len(pdfThrees)):
                for l in range(len(pdfFours)):
                    ind = i+j+k+l
                    prob = pdfOnes[i]*pdfTwos[j]*pdfThrees[k]*pdfFours[l]
                    pdfFinal[ind]+=prob
    return pdfFinal
                    
    

class LandUnitList(object):
    def __init__(self, infantry=0, artillery=0, tanks=0, planes=0, bombers=0, offense=True):
        self.infantry = infantry
        self.artillery = artillery
        self.tanks = tanks
        self.planes=planes
        self.bombers=bombers
        self.offense=offense
        self.update_total()

    def __eq__(self, other):
        try:
            return (self.infantry==other.infantry and \
                    self.artillery==other.artillery and \
                    self.tanks==other.tanks and \
                    self.planes==other.planes and \
                    self.bombers==other.bombers)
        except:
            raise(Exception('Failed equality check on LandUnitList. Possible Type error?'))

    def __str__(self):
        outString = 'Unit List:\n'
        for item in self.__dict__:
            outString+='\t'+item+':\t'+str(self.__dict__[item])+'\n'
        return outString
        
        
    def update_total(self):
        self.totalUnits = self.infantry+self.artillery+self.tanks\
                          +self.planes+self.bombers
        return self.totalUnits

    def get_median_hits(self):
        pdf = self.get_pdf()
        cdf = 0
        for i in range(len(pdf)):
            cdf+=pdf[i]
            if cdf>=0.5:
                return i
        return -1
    
    def get_pdf(self):
        if self.offense:
            return self.get_pdf_offense()
        else:
            return self.get_pdf_defense()
        
    def get_pdf_offense(self):
        nOnes = max(self.infantry-self.artillery,0)
        nTwos = 2*self.artillery - max(self.artillery-self.infantry,0)
        nThrees = self.tanks+self.planes
        nFours = self.bombers

        return generalized_battle_pdf(nOnes, nTwos, nThrees, nFours)

    def get_pdf_defense(self):
        nOnes = self.bombers
        nTwos = self.artillery+self.tanks+self.infantry
        nThrees = 0
        nFours = self.planes

        return generalized_battle_pdf(nOnes, nTwos, nThrees, nFours)

    def allocate_hits(self, nHits, method='normal'):
        tmpUnits = [self.infantry, self.artillery, self.tanks, self.planes, self.bombers]
        if method=='normal':
            if nHits>self.totalUnits:
                return LandUnitList(offense=self.offense)
            i = 0
            while(nHits>0):
                if nHits>=tmpUnits[i]:
                    nHits -= tmpUnits[i]
                    tmpUnits[i]=0
                else:
                    tmpUnits[i] -= nHits
                    nHits=0
            return LandUnitList(infantry=tmpUnits[0], artillery=tmpUnits[1],
                                tanks=tmpUnits[2], planes=tmpUnits[3], bombers=tmpUnits[4],
                                offense=self.offense)
        else:
            raise(Warning('Hit allocation method not defined.'))

        
class SeaUnitList(object):
    def __init__(self, transports=0, subs=0, bships=0, planes=0, bombers=0):
        self.transports = transports
        self.subs = subs
        self.bships = bships
        self.planes = planes
        self.bombers = bombers


class AANode(object):
    def __init__(self, attackingUnits, defendingUnits, landBattle=True):
        self.parentList = []
        self.probability = None
        self.attackingUnits = attackingUnits
        self.defendingUnits = defendingUnits
        self.landBattle = landBattle
        self.update_total()

    def update_total(self):
        self.totalUnits = self.attackingUnits.update_total() + \
                            self.defendingUnits.update_total()
        return self.totalUnits

    def get_probability(self, force=False):
        if not self.probability==None and force==False:
            return self.probability
        
        if self.parentList==[]:
            #print('No parents:')
            #print([self.attackingUnits.infantry, self.defendingUnits.infantry])
            self.probability=1
        else:
            p=0
            for parent, pcgp in self.parentList:
                p+= parent.get_probability()*pcgp
            self.probability=p
        return self.probability
    
    def generate_children(self):
        #FIXME
        pdfAtt = self.attackingUnits.get_pdf()
        pdfDef = self.defendingUnits.get_pdf()
        #must normalize all probabilities to the domain where there exists at least
        #one hit
        normFactor = 1-pdfAtt[0]*pdfDef[0]
        childAttUnits = [[self.attackingUnits.allocate_hits(i), prob]
                         for i, prob in enumerate(pdfDef)]
        childDefUnits = [[self.defendingUnits.allocate_hits(i), prob]
                         for i, prob in enumerate(pdfAtt)]
        #FIXME: remove duplicates from overkill
        
        childList = [[AANode(attUnits, defUnits), aP*dP/normFactor]
                     for attUnits, aP in childAttUnits
                     for defUnits, dP in childDefUnits
                     if not(attUnits==self.attackingUnits and
                            defUnits==self.defendingUnits)]
        return childList

class AANetwork(object):
    def __init__(self):
        self.frontier = []
        self.results = []

    def simulate(self):
        while self.frontier!=[]:
            currentNode = self.frontier.pop(0)
            currentNode.get_probability()
            
            #handle results case
            if (currentNode.attackingUnits == emptyLandUnitList or
                currentNode.defendingUnits == emptyLandUnitList):
                self.results.append(currentNode)
                continue
            
            childList = currentNode.generate_children()
            for child, pcgp in childList:
                #search frontier for matching child. If found, add currentNode to parentList
                #if not found, insert in frontier according to size.
                #
                ###### FIXME ######
                #This search process should be heavily optimized. Frontier
                #is ordered, so we can find where the node should exist,
                #if it does, then insert there. As it is, we're doing
                #two linear searches. Lets go from O(mn) to O(log(mn)) in this
                #section...
                flag = False
                for node in self.frontier:
                    if child.attackingUnits==node.attackingUnits and\
                       child.defendingUnits==node.defendingUnits:
                        flag = True
                        node.parentList.append([currentNode, pcgp])
                        break
                if not flag:
                    child.parentList.append([currentNode, pcgp])
                    numUnits = currentNode.update_total()
                    if self.frontier==[]:
                        self.frontier.append(child)
                    else:
                        for i in range(len(self.frontier)):
                            if (numUnits<=self.frontier[i].update_total() or
                                i==len(self.frontier)-1):
                                self.frontier.insert(i+1,child)
        print('\ncompleted simulation')
        return self.results
                    
                        
emptyLandUnitList = LandUnitList()
        
    

if __name__=='__main__':
    testAttack = LandUnitList(infantry=30)
    testDefense = LandUnitList(infantry=30, offense=False)
    
    rootNode = AANode(testAttack, testDefense)
    
    testNetwork = AANetwork()
    testNetwork.frontier.append(rootNode)

    startTime = time.clock()
    results = testNetwork.simulate()
    endTime = time.clock()
    print('completed in %f seconds' %(endTime-startTime))
    print('\n\n')
    checkSum = 0
    for node in results:
        #print(str(node.attackingUnits)+str(node.defendingUnits)+
        #      "Probability: "+str(node.get_probability()))
        #print('\n')
        print([node.attackingUnits.infantry, node.defendingUnits.infantry,
               node.get_probability()])
        checkSum+=node.get_probability()
    
    print(checkSum)
