# Carlo-Crescente
Cross Section JPSI
//
//  Completiamo quindi con i grafici dell'accettanza e dell'efficienza e con la stima della sezione d'urto differenziale
//

void cross_section(){
    // definiamo i limiti per gli impulsi trasversi delle particelle per  i plot che tracceremo
    double PtMin = 0.;
    double PtMax = 16.;
    // definiamo il numero di bin che vogliamo tracciare nei vari grafici delle popolazioni
    int Nbin = 50;
    
double tagli_pt[51] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.15, 1.3, 1.5, 1.7, 1.95, 2.1, 2.25, 2.5, 2.7, 2.95, 3.1, 3.3, 3.5615, 3.8, 4.0, 4.2, 4.3, 4.54803, 4.63, 4.7, 4.9,   5.15231, 5.358, 5.64166, 5.8, 6.1, 6.35287, 6.78, 7.23158, 7.731, 8.185, 8.30151, 8.54015, 8.76872, 9.22976, 9.72248, 10.325, 11.0662, 11.9326, 12.438, 13.1252, 13.842, 14.3367, 15.2, 16.0};
    
// occupiamoci adesso del calcolo di accettanza ed efficienza, definendo anzitutto gli istogrammi in cui inseriremo i dati di interesse
    TH1D *hg = new TH1D("hg", "Gen JPsi Eta vs Pt", Nbin, PtMin, PtMax);
    TH1D *hrec = new TH1D("hrec", "Rec JPsi Eta vs Pt", Nbin, PtMin, PtMax);
    TH1D *ha = new TH1D("ha", "All good JPsi Eta vs Pt", Nbin, PtMin, PtMax);
    TH1D *hr = new TH1D("hr", "Good taken JPsi Eta vs Pt", Nbin, PtMin, PtMax);
    TH1D *hacc = new TH1D("h_acc", "Accettanza", Nbin, PtMin, PtMax);
    TH1D *heff = new TH1D("h_eff", "Efficienza", Nbin, PtMin, PtMax);
    
// apriamo il file ed accediamo ai dati MC, per l'accettanza e l'efficienza
    TFile *file = new TFile("Jpsi_MC.root");
    TTree *tree = (TTree*) file->Get("data");
    TLorentzVector *gmu1= new TLorentzVector;
    TLorentzVector *gmu2= new TLorentzVector;
    TLorentzVector *gJPsi = new TLorentzVector;
    TLorentzVector *rmu1= new TLorentzVector;
    TLorentzVector *rmu2= new TLorentzVector;
    TLorentzVector *rJPsi = new TLorentzVector;
    
    int Jq;
    int JPsiType;
    int MCType;
    int trig;
   
    tree->SetBranchAddress("muPosP_Gen", &gmu1);
    tree->SetBranchAddress("muNegP_Gen", &gmu2);
    tree->SetBranchAddress("JpsiP_Gen", &gJPsi);

    tree->SetBranchAddress("muPosP", &rmu1);
    tree->SetBranchAddress("muNegP", &rmu2);
    tree->SetBranchAddress("JpsiP", &rJPsi);
    
    tree->SetBranchAddress("JpsiCharge", &Jq);
    tree->SetBranchAddress("JpsiType", &JPsiType);
    tree->SetBranchAddress("MCType", &MCType);
    tree->SetBranchAddress("HLT_DoubleMu0", &trig);
    
    // leggiamo il numero di eventi
    int n = tree->GetEntries();
    cout << "Eventi totali" << n << endl;
    for(int ev = 0; ev < n; ev++) {
        tree->GetEntry(ev);
        
        // leggiamo i dati di interesse per il calcolo dell'accettanza
        double gEta1 = gmu1->Eta();
        double gEta2 = gmu2->Eta();
        double gPt1 = gmu1->Pt();
        double gPt2 = gmu2->Pt();
        double gPtJPsi = gJPsi->Pt();
        double rPtJPsi = rJPsi->Pt();

        // definiamo le variabili booleane che definiscono le condizioni su pseudorapidità ed impulso
        bool gtaglio1 = (fabs(gEta1) < 0.8 && gPt1 > 1.7) || ((fabs(gEta1) > 0.8 && fabs(gEta1) < 1.7) && gPt1 > 0.9) || (fabs(gEta1) > 1.7);
        bool gtaglio2 = (fabs(gEta2) < 0.8 && gPt2 > 1.7) || ((fabs(gEta2) > 0.8 && fabs(gEta2) < 1.7) && gPt2 > 0.9) || (fabs(gEta2) > 1.7);
        hg->Fill(gPtJPsi);
        if(gtaglio1 && gtaglio2){
            if(rJPsi->M() > 3 && rJPsi->M() < 3.2 && Jq == 0 && trig == 1){
                // filliamo gli eventi che soddisfano le condizioni di produzione della JPsi richieste ( = eventi misurati = tutti eventi buoni)
                ha->Fill(rPtJPsi);
                if(MCType == 0 && (JPsiType == 0 || JPsiType == 1 || JPsiType == 2)){
                    // filliamo gli eventi che soddisfanno il tipo di JPsi richiesta ( = eventi buoni)
                    hr->Fill(rPtJPsi);
                }
            }
        }
    }
    hacc->Divide(ha, hg);
    heff->Divide(hr, ha);
      
    TCanvas *c1 = new TCanvas ("c1");
    hacc->Draw();
    hacc->GetXaxis()->SetTitle("p_{T} [ GeV/c ]");
    hacc->GetYaxis()->SetTitle("Acceptance A");
      
    TCanvas *c2 = new TCanvas ("c2");
    heff->Draw();
    heff->GetXaxis()->SetTitle("p_{T} [ GeV/c ]");
    heff->GetYaxis()->SetTitle("Efficiency #varepsilon");
    
    // occupiamoci adesso del calcolo della sezione d'urto
    TFile *file1 = new TFile("Jpsi_DATA.root");
    TTree *tree1 = (TTree*) file1->Get("data");
        
    //  definiamo gli istogrammi per gli impulsi della jpsi e per la sezione d'urto
    TH1D *histoJPsi = new TH1D("J/#Psi p_{T}", "Numero di eventi per bin di p_{T}", Nbin, PtMin, PtMax);
    TH1D *hsigma = new TH1D("h_sigma", "Sezione d'urto differenziale in p_{T}", Nbin, PtMin, PtMax);
    TH1D *hprod = new TH1D("hprod", "hprod", Nbin, PtMin, PtMax);
        
    TLorentzVector *mu1  = new TLorentzVector();
    TLorentzVector *mu2  = new TLorentzVector();
    TLorentzVector *JPsi  = new TLorentzVector();
    int trig1, trig2, Jq1;
    int JPsiType1;

    tree1->SetBranchAddress("HLT_DoubleMu0", &trig1);
    tree1->SetBranchAddress("HLT_DoubleMu0_Quarkonium_v1", &trig2);
    tree1->SetBranchAddress("JpsiCharge", &Jq1);
    tree1->SetBranchAddress("JpsiType", &JPsiType1);
    tree1->SetBranchAddress("muPosP", &mu1);
    tree1->SetBranchAddress("muNegP", &mu2);
    tree1->SetBranchAddress("JpsiP", &JPsi);
    
    int n_tot = tree1->GetEntries();
                
    for (int ev = 0; ev < n_tot; ++ev) {
        tree1->GetEntry(ev);
        double Eta1 = mu1->Eta();
        double Eta2 = mu2->Eta();
        double Pt1 = mu1->Pt();
        double Pt2 = mu2->Pt();
        double PtJPsi = JPsi->Pt();
        bool taglio1 = (fabs(Eta1) < 0.8 && Pt1 > 1.7) || ((fabs(Eta1) > 0.8 && fabs(Eta1) < 1.7) && Pt1 > 0.9) || (fabs(Eta1) > 1.7);
        bool taglio2 = (fabs(Eta2) < 0.8 && Pt2 > 1.7) || ((fabs(Eta2) > 0.8 && fabs(Eta2) < 1.7) && Pt2 > 0.9) || (fabs(Eta2) > 1.7);
        if (JPsi->M() > 3 && JPsi->M() < 3.2  && taglio1 && taglio2 && trig1 == 1 && (JPsiType == 0 || JPsiType == 1 || JPsiType == 2)){
            histoJPsi->Fill(PtJPsi);
        }
    }
        
    TCanvas *c3 = new TCanvas("c3");
    histoJPsi->Draw();
    histoJPsi->GetXaxis()->SetTitle("p_{T} [ GeV/c ]");
    histoJPsi->GetYaxis()->SetTitle("Eventi");
    
    double L = 40 * 1000; //nb^-1
    hprod->Multiply(hacc, heff);
    hsigma->Divide(histoJPsi, hprod);
    hsigma->Scale(1 / L, "width");

    TCanvas *c4 = new TCanvas("c4");
    c4->SetLogy();
    hsigma->Draw();
    hsigma->GetXaxis()->SetTitle("p_{T} [ GeV/c ]");
    hsigma->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}}#upoint BR [ nb/( GeV/c ) ]");
}
