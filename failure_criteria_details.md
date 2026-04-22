# Technical Details: Laminate Failure Criteria

This document provides the mathematical background and implementation strategy for the fracture criteria to be implemented in `CompositCode_main_Teo.m`.

## 1. Maximum Stress Criterion
The Maximum Stress criterion predicts failure when any of the stresses in the principal material directions (local frame 1-2) exceed the corresponding ultimate strength.

### Mathematical Formulation
Failure occurs if:
- **Direction 1 (Fibre):**  
  $\sigma_1 \ge X_T$ (if $\sigma_1 > 0$) or $|\sigma_1| \ge X_C$ (if $\sigma_1 < 0$)
- **Direction 2 (Matrix):**  
  $\sigma_2 \ge Y_T$ (if $\sigma_2 > 0$) or $|\sigma_2| \ge Y_C$ (if $\sigma_2 < 0$)
- **Shear:**  
  $|\tau_{12}| \ge S$

### Failure Index ($IF_{stress}$)
We define a single index where $IF \ge 1$ indicates failure:
$$IF_{stress} = \max \left( \frac{\sigma_1}{X}, \frac{\sigma_2}{Y}, \frac{|\tau_{12}|}{S} \right)$$

---

## 2. Maximum Strain Criterion
Similar to Max Stress, but based on the allowable strains. This is often preferred as it accounts for some Poisson effects that the Max Stress criterion might miss if applied to strains directly.

### Mathematical Formulation
Failure occurs if:
- **Direction 1 (Fibre):**  
  $\epsilon_1 \ge \epsilon_{1T}$ (if $\epsilon_1 > 0$) or $|\epsilon_1| \ge \epsilon_{1C}$ (if $\epsilon_1 < 0$)
- **Direction 2 (Matrix):**  
  $\epsilon_2 \ge \epsilon_{2T}$ (if $\epsilon_2 > 0$) or $|\epsilon_2| \ge \epsilon_{2C}$ (if $\epsilon_2 < 0$)
- **Shear:**  
  $|\gamma_{12}| \ge \gamma_{12S}$

### Failure Index ($IF_{strain}$)
$$IF_{strain} = \max \left( \frac{\epsilon_1}{X_\epsilon}, \frac{\epsilon_2}{Y_\epsilon}, \frac{|\gamma_{12}|}{S_\epsilon} \right)$$

---

## 3. Tsai-Hill Criterion (Quadratic Interaction)
The Tsai-Hill criterion accounts for the interaction between stress components, which is more realistic for composite materials. It is based on the von Mises yield criterion but adapted for orthotropic materials.

### Mathematical Formulation
The Tsai-Hill index is calculated as:
$$\left(\frac{\sigma_1}{X}\right)^2 - \frac{\sigma_1 \sigma_2}{X^2} + \left(\frac{\sigma_2}{Y}\right)^2 + \left(\frac{\tau_{12}}{S}\right)^2 \le 1$$

Where:
- $X$ is $X_T$ if $\sigma_1 > 0$, or $X_C$ if $\sigma_1 < 0$.
- $Y$ is $Y_T$ if $\sigma_2 > 0$, or $Y_C$ if $\sigma_2 < 0$.
- $S$ is the shear strength.

### Implementation Logic
For each lamina, the code will:
1. Determine the signs of $\sigma_1$ and $\sigma_2$.
2. Select the appropriate strengths ($X_T, X_C, Y_T, Y_C$).
3. Compute the quadratic sum.
4. If the sum $\ge 1$, the lamina is considered to have failed.

---

## 4. Implementation Strategy in `CompositCode_main_Teo.m`

### Data Mapping
The `Data` matrix in your file contains the required strengths:
- `Data(6:7, :)`: Sigma1t ($X_T$), Sigma1c ($X_C$)
- `Data(8:9, :)`: Sigma2t ($Y_T$), Sigma2c ($Y_C$)
- `Data(10, :)`: Tau12 ($S$)
- `Data(11:12, :)`: Epsilon1t, Epsilon1c
- `Data(13:14, :)`: Epsilon2t, Epsilon2c
- `Data(15, :)`: Gamma12

### Applied Logic
Since stresses and strains vary linearly through the thickness of each lamina (due to bending moments $M$), we must check **both the bottom and top interfaces** of every lamina.

```matlab
% Conceptual Implementation
for i = 1:num_laminae
    % 1. Get stresses at bot and top for lamina i
    % 2. Calculate IF_max_stress for both, take max
    % 3. Calculate IF_max_strain for both, take max
    % 4. Calculate IF_tsai_hill for both, take max
end
```
