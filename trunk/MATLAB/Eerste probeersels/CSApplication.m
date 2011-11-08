function [] = CSApplication(option)
%CSAPPLICATION Test of CS

t=linspace(0,60,300);
test_function=(sin(t));
percentageVector=[80,50,30,25,10];%
results=[];
error=[];

for d=percentageVector
    if(option==1)
        [res err]=CS(test_function, d);
    elseif(option==2)
        [res err]=CS2(test_function, d);
    else
        disp('Error')
        return;
    end
    results=[results res];
    error=[error err];
end

error
subplot(3,2,1); plot(t,test_function,'k');legend('Original signal')
subplot(3,2,2); plot(t,test_function,'k',t,real(results(:,1)),'b+');legend('CS with 80% samples')
subplot(3,2,3); plot(t,test_function,'k',t,real(results(:,2)),'r+');legend('CS with 50% samples')
subplot(3,2,4); plot(t,test_function,'k',t,real(results(:,3)),'g+');legend('CS with 30% samples')
subplot(3,2,5); plot(t,test_function,'k',t,real(results(:,4)),'y+');legend('CS with 25% samples')
subplot(3,2,6); plot(t,test_function,'k',t,real(results(:,5)),'b+');legend('CS with 10% samples')

end

