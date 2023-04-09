# dt = 0.1
fig1, ax1 = plt.subplots()
ta, ya = anl(0,15,4)
ax1.plot(t,y,'-r',label='Analytical Solution') 
for i in dt:
    t, y = exp(i)
    label = 'Explicit Euler dt=' + str(i) + 's'
    ax1.plot(t,y,'-o',label=label)
ax1.set_title('Explicit Euler Method with Analytical Solution')
ax1.set_xlabel('t')
ax1.set_ylabel('y')
ax1.set_xlim([0,15])
ax1.legend()
plt.savefig('images/0.1_1a.jpg')
plt.show()


# Analytical Solution solution
t_anal = np.linspace(ti,tf,1000)
y_anal = yi * np.exp(-2*t_anal - 0.01/3*t_anal**3)