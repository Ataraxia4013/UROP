for l = 1:range(length(uk(1,:)))
    if all(imag(uk(l,:)) == 0) ==1
        disp('Yes');
        ST = uk(l,:)*Yem;
    end
end