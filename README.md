## Paper: Three Distinct Trajectories of Screen Media Activity in Adolescence: Longitudinal Associations with Accelerated Cortical Thinning and Adverse Behavioral Outcomes

### ***Analysis scripts*** ###

**1. the latent class linear mixed models** (`01_LCGA.R`)

**2. development change scores** (`02_Behaviors_BL4Y.R`)

**3. brain APC** (`03_APC_BL4Y.R`)

**4. normative model** (`04_Normative_BL4Y.R`)

**5. longitduinal environmental factors** (`05_Longitduinal_Environments_2Y4Y.R`)

**6. spatial correlation between t-maps** (`04_Normative_BL4Y_Spin.py`)

**7. bootstrap analysis** (`03_APC_BL4Y_Boot.R`, `04_Normative_BL4Y_Boot.R`)

**8. multi-collinearity** (`06_Multi-collinearity.R`)

### ***Main packages*** ###

**1. R (4.4.3):**
   - **lmerTest** (association analysis)
   = **neuroCombat** (MRI harmonization)
   - **lavvan** (bilateral change score models)

**2. Python (3.12.4):**
   -**neuromap**  (spatial correlation) https://github.com/netneurolab/neuromaps
   -**PCNtookit** (normative model) https://pcntoolkit.readthedocs.io/en/latest/
   -**the pre-trained normative model** https://github.com/predictive-clinical-neuroscience/braincharts