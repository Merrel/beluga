import re
import beluga.bvpsol.algorithms as algorithms
import beluga.optim.Problem
from beluga.optim.problem import *
from beluga.continuation import *

class rashsBackup(object):

    def __init__(self,name='default'):
        self.name = name

        self.independentVariable = []
        self.independentVariableUnits = []

        self.stateVariable = []
        self.stateUnits = []

        self.controlVariable = []
        self.controlUnits = []

        self.costFunctionalExpression = []
        self.costType = []
        self.costUnits = []

        self.dynamicsFunction = []

        self.constraintType = []
        self.constraintExpression = []
        self.constraintUnits = []

        self.constants = []
        self.constantValue = []
        self.constantUnits = []

        self.quantityName = []
        self.quantityExpression = []

        self.solver = algorithms.MultipleShooting(derivative_method='fd',tolerance=1e-4, max_iterations=1000, verbose = True, cached=False, number_arcs=4, max_error=100)

        self.scaleUnits = []
        self.scaleUnitsBy = []

        self.guessMode = []
        self.guessInitialStates = []
        self.guessInitialCoStates = []

        self.mode = 'analytical'

    def independent(self,var,units):
        self.independentVariable = var
        self.independentVariableUnits = units
        return self

    def state(self,var,units):
        self.stateVariable.append(var)
        self.stateUnits.append(units)
        return self

    def control(self,var,units):
        self.controlVariable.append(var)
        self.controlUnits.append(units)
        return self

    def cost(self,type,expr,units):
        self.costType = type
        self.costFunctionalExpression = expr
        self.costUnits = units
        return self

    def dynamics(self,expr):
        self.dynamicsFunction.append(expr)
        return self

    def quantity(self,qty,expr):
        self.quantityName.append(qty)
        self.quantityExpression.append(expr)
        return self

    def constraints(self,type,expr,units):
        self.constraintType.append(type)
        self.constraintExpression.append(expr)
        self.constraintUnits.append(units)
        return self

    def constant(self,var,value,units):
        self.constants.append(var)
        self.constantValue.append(value)
        self.constantUnits.append(units)
        return self

    def scale(self,unit,scaleFactor):
        self.scaleUnits.append(unit)
        self.scaleUnitsBy.append(scaleFactor)
        return self

    def guessSetup(self,guesstype,initialstates,initialcostates):
        self.guessMode = guesstype
        self.guessInitialStates = initialstates
        self.guessInitialCoStates = initialcostates
        return self


