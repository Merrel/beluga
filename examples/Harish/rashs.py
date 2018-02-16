#import re
#import beluga.bvpsol.algorithms as algorithms
#import beluga.optim.Problem
#from beluga.optim.problem import *
#from beluga.continuation import *

class rashs(object):

    def __init__(self,name='default'):
        self.name = name

        self.stateVariable = []
        self.stateUnits = []

        self.pathCostExpression = []
        self.pathCostUnits = []

        self.dynamicsFunction = []

        self.switchingCondition = []

        self.rashsDynamics = []
        self.rashsPathCost = []


    def state(self,var,units):
        self.stateVariable.append(var)
        self.stateUnits.append(units)
        return self

    def pathCost(self,expr,units):
        self.pathCostExpression.append(expr)
        self.pathCostUnits.append(units)
        return self

    def defineDynamicsFunction(self,expr):
        self.dynamicsFunction.append(expr)
        return self

    def switchingConditionSetup(self,expr):
        self.switchingCondition.append(expr)

    def rashsSetup(self):
        self.rashsDynamicsSetup()
        self.rashsPathCostSetup()

        return self


    def rashsDynamicsSetup(self):
        nPhases = len(self.dynamicsFunction[0])
        nStates = len(self.stateVariable)

        for idx1 in range(nStates):
            stateDynamics = ''
            currentStateDynamics = self.dynamicsFunction[idx1]

            for idx2 in range(nPhases):
                sigmoid = '1'
                for idx3 in range(len(self.switchingCondition[idx2])):
                    sigmoidMultiply = '(1/(1+exp((slpRashs*(' + self.switchingCondition[idx2][idx3] + ')))))'
                    sigmoid = sigmoid + '*' + sigmoidMultiply
                dynamicsAttach = '+((' + sigmoid + ')*' + currentStateDynamics[idx2] + ')'
                stateDynamics = stateDynamics + dynamicsAttach
            self.rashsDynamics.append(stateDynamics)
        return self

    def rashsPathCostSetup(self):
        attachCost = ''
        for idx1 in range(len(self.pathCostExpression)):
            sigmoid = '1'
            for idx2 in range(len(self.switchingCondition[idx1])):
                sigmoidMultiply = '(1/(1+exp((slpRashs*(' + self.switchingCondition[idx1][idx2] + ')))))'
                sigmoid = sigmoid + '*' + sigmoidMultiply
            pathCostAttach = '+((' + sigmoid + ')*' + self.pathCostExpression[idx1] + ')'
            attachCost = attachCost + pathCostAttach
        self.rashsPathCost.append(attachCost)
        return self
