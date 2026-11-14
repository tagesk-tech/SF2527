function Q = extrapolate_system_boundaries(Q)
% Continue the nearest interior line to close both artificial boundaries.
Q(:,1) = 2*Q(:,2) - Q(:,3);
Q(:,end) = 2*Q(:,end-1) - Q(:,end-2);
end
