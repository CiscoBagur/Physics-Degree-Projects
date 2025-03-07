# Import the necessary lybraries
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import time


class TwoStateGraph:
    def __init__(self, num_nodes, initial_state_frac, prestige, volatility):
        """
        Initialize a graph with nodes in two possible states.
        
        Args:
        - num_nodes (int): Total number of nodes in the graph
        - initial_state_frac (float): Fraction of initial state being 1
        """
        # Create a random graph with some values for prestige and volatility
        graph = nx.complete_graph(num_nodes)
        self.G = nx.to_numpy_array(graph)
        self.s = prestige
        self.a = volatility
        
        # Initialize node states to the initial fraction(0 or 1)
        self.states = []
        for node in range(num_nodes):
            self.states.append(1 if node/num_nodes < initial_state_frac else 0)
        
        # Store history of state counts
        self.state_history = []
    
    def _count_states(self, num_nodes):
        """
        Count the number of nodes in each state.
        
        Returns:
        - dict: Number of nodes in state 0 and state 1
        """
        total = sum(self.states)
        return total/num_nodes
    
    def update_states(self, update_rule, num_nodes):
        """
        Update node states based on a given update rule.
        
        Args:
        - update_rule (callable): A function that determines a node's next state
          based on its current state and its neighbors' states
        """
        # Create a copy of current states to avoid sequential update bias
        current_states = np.array(self.states.copy())
        # Update each node's state
        for node in range(num_nodes):
            neighbors_index = np.where(self.G[node]==1)[0]
            neighbors = current_states[neighbors_index]
            self.states[node] = update_rule(current_states[node], 
                                            neighbors, prestige=self.s, volatility=self.a)
        
        # Record state counts
        state_counts = self._count_states(num_nodes)
        self.state_history.append(state_counts)
        
        return state_counts
    
    
    def visualize_state_evolution(self):
        """
        Visualize the evolution of node states over time.
        """
        # state_0_history = [counts[0] for counts in self.state_history]
        state_1_history = [counts for counts in self.state_history]
        
        plt.figure(figsize=(10, 6))
        # plt.plot(state_0_history, label='State 0')
        plt.plot(state_1_history, label='Speakers of X')
        plt.title('Node State Evolution')
        plt.xlabel('Timestep')
        plt.ylabel('Number of Nodes')
        plt.ylim(0,1)
        plt.legend()
        plt.show()
        
        
# Rule update
def AS_rule_update(current_state, neighbor_states, prestige, volatility):
    """
    Update rule: node changes language according to the proabilities predicted by the AS model.
    """
    if  sum(neighbor_states)==0:
        return current_state
    
    neighbors_sum = sum(neighbor_states)/len(neighbor_states) # Fraction of neighbours that speak X
    w = np.random.random()
    # print(current_state,neighbor_states,neighbors_sum)
    if current_state==0:
        if w < prestige*neighbors_sum**volatility:
            # print('guanya  x')
            return 1
        else:
            return 0
    if current_state==1:
        if w < (1-prestige)*(1-neighbors_sum)**volatility:
            # print('guanya  y')
            return 0
        else:
            # print('no guanya cap')
            return 1
        
def run_simulation(num_nodes, timesteps, initial_state_frac, prestige, volatility):

    # Calculate the start time
    start = time.time()

    # Create graph
    graph = TwoStateGraph(num_nodes,initial_state_frac, prestige, volatility)
    
    # Initial state counts
    print("Initial state counts:", graph._count_states(num_nodes))
    
    # Loop
    for t in range(timesteps):
        graph.update_states(AS_rule_update,num_nodes)
    
    # Calculate the end time and print it
    end = time.time()
    length = end - start
    print("It took", length, "seconds!")
    
    # Visualize results
    graph.visualize_state_evolution()
    

run_simulation(num_nodes=1000, timesteps=100, initial_state_frac=0.5, prestige=0.55, volatility=0.6)
