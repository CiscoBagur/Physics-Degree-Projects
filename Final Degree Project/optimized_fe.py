import numpy as np
import matplotlib.pyplot as plt
import time


class TwoStateGraph:
    def __init__(self, num_nodes, initial_state_frac, prestige, volatility):
        """
        Initialize a graph with nodes in two possible states.
        
        Args:
        - num_nodes (int): Total number of nodes in the graph
        - initial_state_frac (float): Fraction of initial state being 1
        - prestige (float): Prestige parameter for the AS model
        - volatility (float): Volatility parameter for the AS model
        """
        # Create a complete graph as a NumPy adjacency matrix directly
        self.num_nodes = num_nodes
        self.G = np.ones((num_nodes, num_nodes)) - np.eye(num_nodes)
        self.s = prestige
        self.a = volatility
        
        # Initialize node states more efficiently using NumPy
        self.states = np.zeros(num_nodes)
        initial_ones = int(num_nodes * initial_state_frac)
        self.states[:initial_ones] = 1
        
        # Store history of state counts
        self.state_history = []
        
        # Pre-compute neighbor indices for each node
        self.neighbor_indices = [np.where(self.G[node] == 1)[0] for node in range(num_nodes)]
        
        # Record initial state
        self.state_history.append(self._count_states())
    
    def _count_states(self):
        """
        Count the fraction of nodes in state 1.
        
        Returns:
        - float: Fraction of nodes in state 1
        """
        return np.mean(self.states)
    
    def update_states(self, update_rule):
        """
        Update node states based on a given update rule.
        
        Args:
        - update_rule (callable): A function that determines a node's next state
          based on its current state and its neighbors' states
        """
        # Create a copy of current states to avoid sequential update bias
        current_states = self.states.copy()
        
        # Vectorized random values (generate all at once)
        random_values = np.random.random(self.num_nodes)
        
        # Update each node's state
        for node in range(self.num_nodes):
            neighbors = current_states[self.neighbor_indices[node]]
            self.states[node] = update_rule(current_states[node], 
                                          neighbors, 
                                          self.s, 
                                          self.a, 
                                          random_values[node])
        
        # Record state counts
        state_count = self._count_states()
        self.state_history.append(state_count)
        
        return state_count
    
    def visualize_state_evolution(self):
        """
        Visualize the evolution of node states over time.
        """
        plt.figure(figsize=(10, 6))
        plt.plot(self.state_history, label='Speakers of X')
        plt.title('Node State Evolution')
        plt.xlabel('Timestep')
        plt.ylabel('Fraction of Nodes')
        plt.ylim(0, 1)
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.show()
        
# Optimized rule update
def AS_rule_update(current_state, neighbor_states, prestige, volatility, random_value):
    """
    Update rule: node changes language according to the probabilities predicted by the AS model.
    
    Args:
    - current_state: Current state of the node (0 or 1)
    - neighbor_states: States of neighboring nodes
    - prestige: Prestige parameter
    - volatility: Volatility parameter
    - random_value: Pre-generated random value
    """
    if len(neighbor_states) == 0:
        return current_state
    
    neighbors_sum = np.mean(neighbor_states)  # Fraction of neighbours that speak X
    
    if current_state == 0:
        if random_value < prestige * neighbors_sum**volatility:
            return 1
        else:
            return 0
    else:  # current_state == 1
        if random_value < (1-prestige) * (1-neighbors_sum)**volatility:
            return 0
        else:
            return 1
        
def run_simulation(num_nodes=1000, timesteps=100, initial_state_frac=0.5, prestige=0.55, volatility=0.6, verbose=True):
    """
    Run the AS model simulation.
    
    Args:
    - num_nodes: Number of nodes in the graph
    - timesteps: Number of simulation steps
    - initial_state_frac: Initial fraction of nodes in state 1
    - prestige: Prestige parameter for the AS model
    - volatility: Volatility parameter for the AS model
    - verbose: Whether to print progress information
    """
    # Calculate the start time
    start = time.time()

    # Create graph
    graph = TwoStateGraph(num_nodes, initial_state_frac, prestige, volatility)
    
    # Initial state counts
    if verbose:
        print(f"Initial state counts: {graph.state_history[0]:.3f}")
    
    # Progress tracking
    update_interval = max(1, timesteps // 10)
    
    # Loop
    for t in range(timesteps):
        current_count = graph.update_states(AS_rule_update)
        if verbose and (t+1) % update_interval == 0:
            print(f"Timestep {t+1}/{timesteps}, State 1 fraction: {current_count:.3f}")
    
    # Calculate the end time and print it
    end = time.time()
    runtime = end - start
    if verbose:
        print(f"Simulation completed in {runtime:.2f} seconds!")
    
    # Visualize results
    graph.visualize_state_evolution()
    
    return graph

run_simulation(num_nodes=1000, timesteps=100, initial_state_frac=0.6, prestige=0.55, volatility=0.6)
