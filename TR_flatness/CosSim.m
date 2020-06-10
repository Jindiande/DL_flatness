function sim=CosSim(v1,v2)% cos value of two vector's angle
    sim=abs((transpose(v1)*v2))/(norm(v1,2)*norm(v2,2))
end