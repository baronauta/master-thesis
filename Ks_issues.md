# **$\epsilon=0,\,\Delta=0.2,\,\alpha=0.3,\,T=0.04$**

## Threshold `myMap`

A `myMap` applico `threshold=1e-5`. Notare come in `myMap[392]` compaia **elemento $(2,3)$ non nullo**! Da questo si generano problemi che portano ad avere $Ks$ con termini diagonali non nulli (quando invece dovrebbe essere proporzionale a $\sigma_x$).

![[Screenshot 2025-03-24 114601.png]]

Notare come con `threshold=1e-5` la presenza di termine elemento $(2,3)$ in `myMap` ha come effetto che $Re(Ks_{00}) \neq 0$; inoltre si ha $Im(Ks_{01}) \neq 0$.

![[Screenshot 2025-03-25 081051.png]]
![[Screenshot 2025-03-25 081059 1.png]]

---
Alzando il valore di soglia a `threshold=1e-4` riesco ad eliminare la presenza del termine $(2,3)$ (ma solo in questo particolare set di parametri): segue che  $Re(Ks_{00})=0$ e $Im(Ks_{01})=0$ (quasi).

![[Screenshot 2025-03-25 081129 1.png]]
![[Screenshot 2025-03-25 081142.png]]

---

Se scelgo valore di soglia troppo alto, e.g. `threshold=1e-2`, poi non riesco a calcolare `myGen`. Prendendo `threshold=1e-3` riesco a procedere nei calcoli ma questo non mi evita il comportamento indesiderato in `myMap`: ho comunque il termine $(2,3)$ non nullo.

## Threshold `meas`

**Threshold su misurazioni** è ininfluente!

## Threshold `eigvals` e `eigvecs`

Threshold su `eigvals` e `eigenvecs` necessario per visualizzazione appropriata dei risultati, altrimenti ho termini per $Ks$ sostanzialmente nulli, e.g $1e-17$, che creano difficoltà grafiche con `Plots`. Questo threshold è sostanzialmente per **esigenze grafiche**. Quello per `eigvals` è inutile dal momento che non vengono interessati dal threshold: autovalori sono nell'ordine di $1e-1$. Threshold per `eigvecs` serve, ma si può scegliere anche indipendentemente da quello per `myMap`.


# **$\epsilon=0,\,\Delta=0.2,\,\alpha=0.1,\,T=0.04$

Cambiando set di parametri la scelta di `threshold` è nuovamente problematica: non solo si ha  $Re(Ks_{00}) \neq 0$, ma si ha anche la presenza di picchi anomali in $Re(Ks_{01})$.

---
Con `threshold = 1e-4` ho problemi per il termine $Re(Ks_{01})$ per tempi più lunghi. In particolare compaiono dei "picchetti" anomali. Voglio capire cosa li genera. 

![[Screenshot 2025-03-24 122054.png]]

Studio il primo **picco anomalo** che si osserva, mi sto riferendo a quello attorno a $t=12.5$.
Per quanto riguarda `myMap` noto che i termini corrispondenti al picco sono quelli dove il **termine $(4,3)$ diventa nullo**!

![[Screenshot 2025-03-25 091403 1.png]]

Situazione analoga anche per il picco successivo, i.d. $t=16$, dove ad essere nullo è il termine $(3,3)$.

![[Screenshot 2025-03-25 092335.png]]

Analogo comportamento per il picco intorno a $t=21$ dove ora ad essere nullo è il termine $(4,4)$.

---

Con `threshold = 1e-5` ho meno picchi per $Re(Ks_{01})$ rispetto al caso precedente. Perchè si generano?
Stessa situazione di prima: uno dei termini nel blocco $2 \times 2$ in basso a destra diventa nullo.

![[Screenshot 2025-03-24 121732 1.png]]

La differenza con il caso precedente è che con questo threshold più alto sopravvivono più termini in questo blocco (che hanno valori attorno al valore di soglia) e dunque non venendo forzatamente messi a zero, non si ha la generazione del picco. D'altra parte l'imposizione del threshold è necessaria perché altrimenti ho elementi non nulli nel blocco in alto a destra che producono $Re(Ks_{00}) \neq 0$ come si vede nel caso di `threshold=1e-6`.

---

Con `threshold=1e-6`, il termine $Re(Ks_{01})$ si regolarizza tuttavia si "sporca" la sua componente immaginare e $Ks_{00}$ (che nei precedenti casi era esattamente nullo, come dovrebbe) ricompare in modo problematico.

![[Screenshot 2025-03-24 123004.png]]![[Screenshot 2025-03-24 122810 1.png]]

---

Impongo blocco $2\times2$ in alto a dx e in basso a sx nulli e non impongo alcun threshold ulteriore su `myMap`:
- $Re(Ks_{00})=0$ come deve
- Non ho picchi anomali in $Re(Ks_{01})$
- **MA** continuo ad avere comportamento oscillante $Re(Ks_{01})$ che non è presente nei plot del paper *Gatto et al.*
![[Screenshot 2025-03-25 100242.png]]![[Screenshot 2025-03-25 100259.png]]

Imponendo poi il threshold a `myMap`, se il threshold interessa i valori negli altri due blocchi allora si generano i picchetti anomali in $Re(Ks_{01})$.